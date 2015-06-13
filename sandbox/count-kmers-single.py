#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) University of California, Davis, 2015. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=missing-docstring,invalid-name
"""
Produce k-mer counts for all the k-mers in the given sequence file,
using the given counting table.

% python sandbox/count-kmers-single.py <fasta/fastq>

Use '-h' for parameter help.
"""
from __future__ import print_function

import sys
import khmer
import argparse
import screed
import csv
from khmer.khmer_args import (build_counting_args, report_on_config, info,
                              add_threading_args)
from khmer.kfile import (check_input_files, check_space,
                         check_space_for_hashtable)
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

    print ('making k-mer counting table', file=sys.stderr)
    counting_hash = khmer.new_counting_hash(args.ksize, args.min_tablesize,
                                            args.n_tables)
    # @CTB counting_hash.set_use_bigcount(args.bigcount)

    kmer_size = counting_hash.ksize()
    hashsizes = counting_hash.hashsizes()
    tracking = khmer._Hashbits(  # pylint: disable=protected-access
        kmer_size, hashsizes)

    print ('kmer_size: %s' % counting_hash.ksize(), file=sys.stderr)
    print ('k-mer counting table sizes: %s' % (counting_hash.hashsizes(),),
           file=sys.stderr)

    if args.output_file is None:
        args.output_file = sys.stdout
    writer = csv.writer(args.output_file)

    # start loading
    rparser = khmer.ReadParser(args.input_sequence_filename)
    threads = []
    print ('consuming input, round 1 -- %s' % (args.input_sequence_filename),
           file=sys.stderr)
    for _ in xrange(args.threads):
        thread = \
            threading.Thread(
                target=counting_hash.consume_fasta_with_reads_parser,
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
                writer.writerow([kmer, str(counting_hash.get(kmer))])

    print ('Total number of unique k-mers: {0}'.format(
        counting_hash.n_unique_kmers()), file=sys.stderr)


if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
