#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=missing-docstring,invalid-name
"""
Count the median/avg k-mer abundance for each sequence in the input file,
based on the k-mer counts in the given k-mer counting table.  Can be used to
estimate expression levels (mRNAseq) or coverage (genomic/metagenomic).

% scripts/count-median.py <htname> <input seqs> <output counts>

Use '-h' for parameter help.

The output file contains sequence id, median, average, stddev, and seq length.

NOTE: All 'N's in the input sequences are converted to 'G's.
"""
import screed
import argparse
import khmer
from khmer.file import check_file_status, check_space
import textwrap


def get_parser():
    epilog = """
    Count the median/avg k-mer abundance for each sequence in the input file,
    based on the k-mer counts in the given k-mer counting table.  Can be used
    to estimate expression levels (mRNAseq) or coverage (genomic/metagenomic).

    The output file contains sequence id, median, average, stddev, and seq
    length.

    NOTE: All 'N's in the input sequences are converted to 'G's.
    """
    parser = argparse.ArgumentParser(
        description='Count k-mers summary stats for sequences',
        epilog=textwrap.dedent(epilog))

    parser.add_argument('ctfile', help='input k-mer count table filename')
    parser.add_argument('input', help='input FAST[AQ] sequence filename')
    parser.add_argument('output', help='output summary filename')
    parser.add_argument('--version', action='version', version='%(prog)s '
                        + khmer.__version__)
    return parser


def main():
    args = get_parser().parse_args()

    htfile = args.ctfile
    input_filename = args.input
    output_filename = args.output

    infiles = [htfile, input_filename]
    for infile in infiles:
        check_file_status(infile)

    check_space(infiles)

    print 'loading k-mer counting table from', htfile
    htable = khmer.load_counting_hash(htfile)
    ksize = htable.ksize()

    print 'writing to', output_filename
    output = open(output_filename, 'w')

    for record in screed.open(input_filename):
        seq = record.sequence.upper()
        if 'N' in seq:
            seq = seq.replace('N', 'G')

        if ksize <= len(seq):
            medn, ave, stdev = htable.get_median_count(seq)
            print >> output, record.name, medn, ave, stdev, len(seq)

if __name__ == '__main__':
    main()
