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

% python sandbox/count-kmers-single.py <fasta/fastq>

Use '-h' for parameter help.
"""

import sys
import khmer
import argparse
import screed
import csv
from khmer.khmer_args import (build_counting_args, report_on_config, info,
                              add_threading_args)
from khmer.kfile import (check_input_files, check_space,
                         check_space_for_graph)
import threading


def get_parser():
    parser = build_counting_args(
        descr="Output abundances of the k-mers in the sequence file.")
    add_threading_args(parser)

    parser.add_argument('input_sequence_filename', help='The input'
                        ' FAST[AQ] sequence file.')

    parser.add_argument('-o', '--out', metavar="output_file",
                        dest='output_file',
                        type=argparse.FileType('w'),
                        default=None, help='output counts to this file')

    return parser


def main():
    info('count-kmers-single.py', ['counting'])
    args = get_parser().parse_args()

    check_input_files(args.input_sequence_filename, False)

    print ('making k-mer countgraph', file=sys.stderr)
    countgraph = khmer.Countgraph(args.ksize, args.max_tablesize,
                                            args.n_tables)
    # @CTB countgraph.set_use_bigcount(args.bigcount)

    kmer_size = countgraph.ksize()
    hashsizes = countgraph.hashsizes()
    tracking = khmer.Nodegraph(  # pylint: disable=protected-access
        kmer_size, 1, 1, primes=hashsizes)

    print ('kmer_size: %s' % countgraph.ksize(), file=sys.stderr)
    print ('k-mer countgraph sizes: %s' % (countgraph.hashsizes(),),
           file=sys.stderr)

    if args.output_file is None:
        args.output_file = sys.stdout
    writer = csv.writer(args.output_file)

    # start loading
    rparser = khmer.ReadParser(args.input_sequence_filename)
    threads = []
    print ('consuming input, round 1 -- %s' % (args.input_sequence_filename),
           file=sys.stderr)
    for _ in range(args.threads):
        thread = \
            threading.Thread(
                target=countgraph.consume_seqfile,
                args=(rparser, )
            )
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()

    for record in screed.open(args.input_sequence_filename):
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
