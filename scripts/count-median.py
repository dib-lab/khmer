#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2012-2015, Michigan State University.
# Copyright (C) 2015, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org
# pylint: disable=missing-docstring,invalid-name
"""
Count the median/avg k-mer abundance for each sequence in the input file.

The abundance is based on the k-mer counts in the given k-mer countgraph.
Can be used to estimate expression levels (mRNAseq) or coverage
(genomic/metagenomic).

% scripts/count-median.py <countgraph> <input seqs> <output counts>

Use '-h' for parameter help.

The output file contains sequence id, median, average, stddev, and seq length.

NOTE: All 'N's in the input sequences are converted to 'A's.
"""
import argparse
import screed
import sys
import csv
import textwrap

from khmer import __version__, Countgraph
from khmer.kfile import check_input_files, check_space
from khmer.khmer_args import (sanitize_help, KhmerArgumentParser)


def get_parser():
    epilog = """\
    Count the median/avg k-mer abundance for each sequence in the input file,
    based on the k-mer counts in the given k-mer countgraph.  Can be used
    to estimate expression levels (mRNAseq) or coverage (genomic/metagenomic).

    The output file contains sequence id, median, average, stddev, and
    seq length, in comma-separated value (CSV) format.

    Example::

        load-into-counting.py counts tests/test-data/test-reads.fq.gz
        count-median.py counts tests/test-data/test-reads.fq.gz medians.txt

    NOTE: All 'N's in the input sequences are converted to 'A's.
    """
    parser = KhmerArgumentParser(
        description='Count k-mers summary stats for sequences',
        epilog=textwrap.dedent(epilog))

    parser.add_argument('countgraph', metavar='input_count_graph_filename',
                        help='input k-mer countgraph filename')
    parser.add_argument('input', metavar='input_sequence_filename',
                        help='input FAST[AQ] sequence filename')
    parser.add_argument('output', metavar='output_summary_filename',
                        help='output summary filename',
                        type=argparse.FileType('w'))
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Overwrite output file if it exists')
    return parser


def main():
    args = sanitize_help(get_parser()).parse_args()

    htfile = args.countgraph
    input_filename = args.input
    output = args.output

    infiles = [htfile, input_filename]
    for infile in infiles:
        check_input_files(infile, args.force)

    check_space(infiles, args.force)

    print('loading k-mer countgraph from', htfile, file=sys.stderr)
    countgraph = Countgraph.load(htfile)
    ksize = countgraph.ksize()
    print('writing to', output.name, file=sys.stderr)

    output = csv.writer(output)
    # write headers:
    output.writerow(['name', 'median', 'average', 'stddev', 'seqlen'])

    for record in screed.open(input_filename):
        seq = record.sequence.upper()
        if 'N' in seq:
            seq = seq.replace('N', 'A')

        if ksize <= len(seq):
            medn, ave, stdev = countgraph.get_median_count(seq)
            ave, stdev = [round(x, 9) for x in (ave, stdev)]
            output.writerow([record.name, medn, ave, stdev, len(seq)])


if __name__ == '__main__':
    main()
