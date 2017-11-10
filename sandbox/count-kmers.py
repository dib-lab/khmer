#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
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
Produce k-mer counts for all the k-mers in the given sequence file,
using the given countgraph.

% python sandbox/count-kmers.py <ct> <fasta/fastq> [ <fasta/fastq> ... ]

Use '-h' for parameter help.
"""

import sys
import khmer
import argparse
import screed
import csv
from khmer import Countgraph
from khmer.khmer_args import info


def get_parser():
    parser = argparse.ArgumentParser(
        description="Output abundances of the k-mers in "
        "the sequence files using a pre-made k-mer countgraph.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input_count_graph_filename', help='The name of the'
                        ' input k-mer countgraph file.')
    parser.add_argument('input_sequence_filenames', help='The input'
                        ' FAST[AQ] sequence file(s).', nargs='+')

    parser.add_argument('-o', '--out', metavar="output_file",
                        dest='output_file',
                        type=argparse.FileType('w'),
                        default=None, help='output counts to this file')

    return parser


def main():
    info('count-kmers.py', ['counting'])
    args = get_parser().parse_args()

    print ('hashtable from', args.input_count_graph_filename,
           file=sys.stderr)
    countgraph = Countgraph.load(
        args.input_count_graph_filename)

    kmer_size = countgraph.ksize()
    hashsizes = countgraph.hashsizes()
    tracking = khmer.Nodegraph(  # pylint: disable=protected-access
        kmer_size, 1, 1, primes=hashsizes)

    if args.output_file is None:
        args.output_file = sys.stdout
    writer = csv.writer(args.output_file)

    for filename in args.input_sequence_filenames:
        for record in screed.open(filename):
            seq = record.sequence.replace('N', 'A')
            for i in range(len(seq) - kmer_size + 1):
                kmer = seq[i:i+kmer_size]
                if not tracking.get(kmer):
                    tracking.count(kmer)
                    writer.writerow([kmer, str(countgraph.get(kmer))])

    print ('Total number of unique k-mers: {0}'.format(
        countgraph.n_unique_kmers()), file=sys.stderr)


if __name__ == '__main__':
    main()

# vim: set filetype=python tabstop=4 softtabstop=4 shiftwidth=4 expandtab:
# vim: set textwidth=79:
