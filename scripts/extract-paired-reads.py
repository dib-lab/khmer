#! /usr/bin/env python
#
# This script is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring
"""
Take a file containing a mixture of interleaved and orphaned reads, and
extract them into separate files (.pe and .se).

% scripts/extract-paired-reads.py <infile>

Reads FASTQ and FASTA input, retains format for output.
"""
import screed
import sys
import os.path
import textwrap
import argparse
import khmer
from khmer.file import check_file_status, check_space


def is_pair(name1, name2):
    if name1.endswith('/1') and name2.endswith('/2'):
        subpart1 = name1.split('/')[0]
        subpart2 = name2.split('/')[0]
        if subpart1 == subpart2:
            assert subpart1
            return True

    return False


def output_pair(read1, read2):
    if hasattr(read1, 'accuracy'):
        return "@%s\n%s\n+\n%s\n@%s\n%s\n+\n%s\n" % \
            (read1.name, read1.sequence, read1.accuracy,
             read2.name, read2.sequence, read2.accuracy)
    else:
        return ">%s\n%s\n>%s\n%s\n" % (read1.name, read1.sequence, read2.name,
                                       read2.sequence)


def output_single(read):
    if hasattr(read, 'accuracy'):
        return "@%s\n%s\n+\n%s\n" % (read.name, read.sequence, read.accuracy)
    else:
        return ">%s\n%s\n" % (read.name, read.sequence)


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
        description=
        'Take a mixture of reads and split into pairs and orphans.',
        epilog=textwrap.dedent(epilog))
    parser.add_argument('infile')
    parser.add_argument('--version', action='version', version='%(prog)s '
                        + khmer.__version__)
    return parser


def main():
    args = get_parser().parse_args()

    check_file_status(args.infile)
    infiles = [args.infile]
    check_space(infiles)

    outfile = os.path.basename(args.infile)
    if len(sys.argv) > 2:
        outfile = sys.argv[2]

    single_fp = open(outfile + '.se', 'w')
    paired_fp = open(outfile + '.pe', 'w')

    print 'reading file "%s"' % args.infile
    print 'outputting interleaved pairs to "%s.pe"' % outfile
    print 'outputting orphans to "%s.se"' % outfile

    last_record = None
    last_name = None

    n_pe = 0
    n_se = 0

    record = None
    for index, record in enumerate(screed.open(sys.argv[1])):
        if index % 100000 == 0 and index > 0:
            print '...', index
        name = record['name'].split()[0]

        if last_record:
            if is_pair(last_name, name):
                paired_fp.write(output_pair(last_record, record))
                name, record = None, None
                n_pe += 1
            else:
                single_fp.write(output_single(last_record))
                n_se += 1

        last_name = name
        last_record = record

    if last_record:
        if is_pair(last_name, name):
            paired_fp.write(output_pair(last_record, record))
            name, record = None, None
            n_pe += 1
        else:
            single_fp.write(output_single(last_record))
            name, record = None, None
            n_se += 1

    if record:
        single_fp.write(output_single(record))
        n_se += 1

    single_fp.close()
    paired_fp.close()

    if n_pe == 0:
        raise Exception("no paired reads!? check file formats...")

    print 'DONE; read %d sequences, %d pairs and %d singletons' % \
          (index + 1, n_pe, n_se)

if __name__ == '__main__':
    main()
