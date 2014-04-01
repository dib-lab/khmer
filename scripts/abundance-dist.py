#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2010-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=missing-docstring,invalid-name
"""
Produce the k-mer abundance distribution for the given file.

% python scripts/abundance-dist.py [ -z -s ] <htname> <data> <histout>

Use '-h' for parameter help.
"""
from __future__ import print_function

import sys
import khmer
import argparse
import os
from khmer.file import check_file_status, check_space
from khmer.khmer_args import info


def get_parser():
    parser = argparse.ArgumentParser(
        description="Calculate abundance distribution of the k-mers in "
        "the sequence file using a pre-made k-mer counting table.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input_counting_table_filename', help='The name of the'
                        ' input k-mer counting table file.')
    parser.add_argument('input_sequence_filename', help='The name of the input'
                        ' FAST[AQ] sequence file.')
    parser.add_argument('output_histogram_filename', help='The columns are: '
                        '(1) k-mer abundance, (2) k-mer count, (3) cumulative '
                        'count, (4) fraction of total distinct k-mers.')
    parser.add_argument('-z', '--no-zero', dest='output_zero', default=True,
                        action='store_false',
                        help='Do not output 0-count bins')
    parser.add_argument('-s', '--squash', dest='squash_output', default=False,
                        action='store_true',
                        help='Overwrite output file if it exists')
    parser.add_argument('--version', action='version', version='%(prog)s '
                        + khmer.__version__)
    return parser


def main():
    info('abundance-dist.py', ['counting'])
    args = get_parser().parse_args()
    infiles = [args.input_counting_table_filename,
               args.input_sequence_filename]
    for infile in infiles:
        check_file_status(infile)

    check_space(infiles)

    print('hashtable from', args.input_counting_table_filename)
    counting_hash = khmer.load_counting_hash(
        args.input_counting_table_filename)

    kmer_size = counting_hash.ksize()
    hashsizes = counting_hash.hashsizes()
    tracking = khmer._new_hashbits(  # pylint: disable=protected-access
        kmer_size, hashsizes)

    print('K:', kmer_size)
    print('HT sizes:', hashsizes)
    print('outputting to', args.output_histogram_filename)

    if os.path.exists(args.output_histogram_filename):
        if not args.squash_output:
            print('ERROR: %s exists; not squashing.' %
                  args.output_histogram_filename,
                  file=sys.stderr)
            sys.exit(1)

        print('** squashing existing file %s' % args.output_histogram_filename)

    print('preparing hist...')
    abundances = counting_hash.abundance_distribution(
        args.input_sequence_filename, tracking)
    total = sum(abundances)

    if 0 == total:
        print("ERROR: abundance distribution is uniformly zero; "
              "nothing to report.", file=sys.stderr)
        print("\tPlease verify that the input files are valid.",
              file=sys.stderr)
        sys.exit(1)
    hash_fp = open(args.output_histogram_filename, 'w')

    sofar = 0
    for _, i in enumerate(abundances):
        if i == 0 and not args.output_zero:
            continue

        sofar += i
        frac = sofar / float(total)

        print(_, i, sofar, round(frac, 3), file=hash_fp)

        if sofar == total:
            break

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
