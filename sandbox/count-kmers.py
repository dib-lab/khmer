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

% python sandbox/count-kmers.py <ct> <fasta/fastq> [ <fasta/fastq> ... ]

Use '-h' for parameter help.
"""
from __future__ import print_function

import sys
import khmer
import argparse
import screed
import csv
from khmer.khmer_args import info

# TODO
# write 'count-kmers-single.py' a la abundance-dist-single

def get_parser():
    parser = argparse.ArgumentParser(
        description="Output abundances of the k-mers in "
        "the sequence file using a pre-made k-mer counting table.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('input_counting_table_filename', help='The name of the'
                        ' input k-mer counting table file.')
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

    print ('hashtable from', args.input_counting_table_filename,
           file=sys.stderr)
    counting_hash = khmer.load_counting_hash(
        args.input_counting_table_filename)
    
    kmer_size = counting_hash.ksize()
    hashsizes = counting_hash.hashsizes()
    tracking = khmer._Hashbits(  # pylint: disable=protected-access
        kmer_size, hashsizes)

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
                    writer.writerow([kmer, str(counting_hash.get(kmer))])
    

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
