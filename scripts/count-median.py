#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# pylint: disable=missing-docstring,invalid-name
"""
Count the median/avg k-mer abundance for each sequence in the input file.

The abundance is based on the k-mer counts in the given k-mer counting
table.  Can be used to estimate expression levels (mRNAseq) or coverage
(genomic/metagenomic).

% scripts/count-median.py <htname> <input seqs> <output counts>

Use '-h' for parameter help.

The output file contains sequence id, median, average, stddev, and seq length.

NOTE: All 'N's in the input sequences are converted to 'A's.
"""
from __future__ import print_function
import screed
import argparse
import sys
import csv
import textwrap

import khmer
from khmer.kfile import check_input_files, check_space
from khmer.khmer_args import info


def get_parser():
    epilog = """
    Count the median/avg k-mer abundance for each sequence in the input file,
    based on the k-mer counts in the given k-mer counting table.  Can be used
    to estimate expression levels (mRNAseq) or coverage (genomic/metagenomic).

    The output file contains sequence id, median, average, stddev, and
    seq length; fields are separated by spaces. For khmer 1.x
    count-median.py will split sequence names at the first space which
    means that some sequence formats (e.g. paired FASTQ in Casava 1.8
    format) will yield uninformative names.  Use :option:`--csv` to
    fix this behavior.

    Example::

       count-median.py counts.ct tests/test-data/test-reads.fq.gz medians.txt

    NOTE: All 'N's in the input sequences are converted to 'A's.
    """
    parser = argparse.ArgumentParser(
        description='Count k-mers summary stats for sequences',
        epilog=textwrap.dedent(epilog))

    parser.add_argument('ctfile', metavar='input_counting_table_filename',
                        help='input k-mer count table filename')
    parser.add_argument('input', metavar='input_sequence_filename',
                        help='input FAST[AQ] sequence filename')
    parser.add_argument('output', metavar='output_summary_filename',
                        help='output summary filename')
    parser.add_argument('--version', action='version', version='%(prog)s ' +
                        khmer.__version__)
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Overwrite output file if it exists')
    parser.add_argument('--csv', default=False, action='store_true',
                        help="Use the CSV format for the histogram."
                        "Includes column headers.")
    return parser


def main():
    info('count-median.py', ['diginorm'])
    args = get_parser().parse_args()

    htfile = args.ctfile
    input_filename = args.input
    output_filename = args.output

    infiles = [htfile, input_filename]
    for infile in infiles:
        check_input_files(infile, args.force)

    check_space(infiles, args.force)

    print('loading k-mer counting table from', htfile, file=sys.stderr)
    htable = khmer.load_counting_hash(htfile)
    ksize = htable.ksize()

    print('writing to', output_filename, file=sys.stderr)
    output = open(output_filename, 'w')

    if args.csv:
        output = csv.writer(output)
        # write headers:
        output.writerow(['name', 'median', 'average', 'stddev', 'seqlen'])

    parse_description = True            # @legacy behavior: split seq headers
    if args.csv:
        parse_description = False       # only enable if we're doing csv out

    for record in screed.open(input_filename,
                              parse_description=parse_description):
        seq = record.sequence.upper()
        if 'N' in seq:
            seq = seq.replace('N', 'A')

        if ksize <= len(seq):
            medn, ave, stdev = htable.get_median_count(seq)
            ave, stdev = [round(x, 9) for x in (ave, stdev)]
            if args.csv:
                output.writerow([record.name, medn, ave, stdev, len(seq)])
            else:
                print(record.name, medn, ave, stdev, len(seq), file=output)

if __name__ == '__main__':
    main()
