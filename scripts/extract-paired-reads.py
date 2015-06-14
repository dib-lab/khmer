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
    Separate properly interleaved reads from orphaned reads.

    Many read-handling programs (assemblers, mappers, etc.) require
    that you give them either perfectly interleaved files, or files
    containing only single reads. This script takes files that were
    originally interleaved but where reads may have been orphaned (via
    error filtering, application of abundance filtering, digital
    normalization in non-paired mode, or partitioning) and separates
    the interleaved reads from the orphaned reads.

    The default output is two files, <input file>.pe and <input
    file>.se, placed in the current directory. The .pe file contains
    interleaved and properly paired sequences, while the .se file
    contains orphan sequences.

    The directory into which the interleaved and orphaned reads are
    output may be specified using :option:`-o`/:option:`--output-dir`.
    This directory will be created if it does not already exist.

    Alternatively, you can specify the filenames directly with
    :option:`-p`/:option:`--output-paired` and
    :option:`-s`/:option:`--output-single`, which will override the
    :option:`-o`/:option:`--output-dir` option.

    Example::

        extract-paired-reads.py tests/test-data/paired.fq
    """
    parser = argparse.ArgumentParser(
        description='Take a mixture of reads and split into pairs and '
        'orphans.', epilog=textwrap.dedent(epilog))
    parser.add_argument('infile')
    parser.add_argument('--version', action='version', version='%(prog)s ' +
                        khmer.__version__)

    parser.add_argument('-o', '--output-dir', metavar="output_directory",
                        dest='output_directory', default='', help='Output '
                        'split reads to specified directory. Creates '
                        'directory if necessary')

    parser.add_argument('-p', '--output-paired', metavar='output_paired',
                        default=None, help='Output paired reads to this '
                        'file', type=argparse.FileType('w'))
    parser.add_argument('-s', '--output-single', metavar='output_single',
                        default=None, help='Output orphaned reads to this '
                        'file', type=argparse.FileType('w'))

    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Overwrite output file if it exists')
    return parser


def main():
    info('extract-paired-reads.py')
    args = get_parser().parse_args()

    infile = args.infile
    check_input_files(infile, args.force)
    check_space([infile], args.force)

    # decide where to put output files - specific directory? or just default?
    if infile == '/dev/stdin':
        if not (args.output_paired and args.output_single):
            print("Accepting input from stdin; output filenames must be "
                  "provided.", file=sys.stderr)
            sys.exit(1)
    elif args.output_directory:
        if not os.path.exists(args.output_directory):
            os.makedirs(args.output_directory)
        out1 = args.output_directory + '/' + os.path.basename(infile) + '.se'
        out2 = args.output_directory + '/' + os.path.basename(infile) + '.pe'
        single_fp = open(out1, 'w')
        paired_fp = open(out2, 'w')
    else:
        out1 = os.path.basename(infile) + '.se'
        out2 = os.path.basename(infile) + '.pe'
        single_fp = open(out1, 'w')
        paired_fp = open(out2, 'w')

    # OVERRIDE default output file locations with -p, -s
    if args.output_paired:
        paired_fp = args.output_paired
        out2 = paired_fp.name
    if args.output_single:
        single_fp = args.output_single
        out1 = single_fp.name

    print('reading file "%s"' % infile, file=sys.stderr)
    print('outputting interleaved pairs to "%s"' % out2, file=sys.stderr)
    print('outputting orphans to "%s"' % out1, file=sys.stderr)

    n_pe = 0
    n_se = 0

    screed_iter = screed.open(infile, parse_description=False)
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

    print('wrote to: %s and %s' % (out2, out1),
          file=sys.stderr)


if __name__ == '__main__':
    main()
