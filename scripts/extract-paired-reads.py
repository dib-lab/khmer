#! /usr/bin/env python
#
# This script is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring
"""
Split up pairs and singletons.

Take a file containing a mixture of interleaved and orphaned reads, and
extract them into separate files (.pe and .se).

% scripts/extract-paired-reads.py <infile>

Reads FASTQ and FASTA input, retains format for output.
"""
from __future__ import print_function
import screed
import sys
import os.path
import textwrap
import argparse
import khmer
from khmer.kfile import check_input_files, check_space
from khmer.khmer_args import info

from khmer.utils import broken_paired_reader, write_record, write_record_pair


def get_parser():
    epilog = """
    The output is two files, <input file>.pe and <input file>.se, placed in the
    current directory. The .pe file contains interleaved and properly paired
    sequences, while the .se file contains orphan sequences.

    Many assemblers (e.g. Velvet) require that you give them either perfectly
    interleaved files, or files containing only single reads. This script takes
    files that were originally interleaved but where reads may have been
    orphaned via error filtering, application of abundance filtering, digital
    normalization in non-paired mode, or partitioning.

    Example::

        extract-paired-reads.py tests/test-data/paired.fq
    """
    parser = argparse.ArgumentParser(
        description='Take a mixture of reads and split into pairs and '
        'orphans.', epilog=textwrap.dedent(epilog))
    parser.add_argument('infile')
    parser.add_argument('--version', action='version', version='%(prog)s ' +
                        khmer.__version__)
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Overwrite output file if it exists')
    return parser


def main():
    info('extract-paired-reads.py')
    args = get_parser().parse_args()

    check_input_files(args.infile, args.force)
    infiles = [args.infile]
    check_space(infiles, args.force)

    outfile = os.path.basename(args.infile)
    if len(sys.argv) > 2:
        outfile = sys.argv[2]

    single_fp = open(outfile + '.se', 'w')
    paired_fp = open(outfile + '.pe', 'w')

    print('reading file "%s"' % args.infile, file=sys.stderr)
    print('outputting interleaved pairs to "%s.pe"' % outfile, file=sys.stderr)
    print('outputting orphans to "%s.se"' % outfile, file=sys.stderr)

    n_pe = 0
    n_se = 0

    screed_iter = screed.open(args.infile, parse_description=False)
    for index, is_pair, read1, read2 in broken_paired_reader(screed_iter):
        if index % 100000 == 0 and index > 0:
            print('...', index, file=sys.stderr)

        if is_pair:
            write_record_pair(read1, read2, paired_fp)
            n_pe += 1
        else:
            write_record(read1, single_fp)
            n_se += 1

    single_fp.close()
    paired_fp.close()

    if n_pe == 0:
        raise Exception("no paired reads!? check file formats...")

    print('DONE; read %d sequences,'
          ' %d pairs and %d singletons' %
          (n_pe * 2 + n_se, n_pe, n_se), file=sys.stderr)

    print('wrote to: ' + outfile + '.se' + ' and ' + outfile + '.pe',
          file=sys.stderr)


if __name__ == '__main__':
    main()
