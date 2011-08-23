#! /usr/bin/env python
"""
Eliminate reads with median k-mer abundance higher than DESIRED_COVERAGE.

Parameters to adjust: DESIRED_COVERAGE, 
"""

import sys, screed, os
import khmer
import random

DESIRED_COVERAGE=20

import argparse

DEFAULT_K=32
DEFAULT_N_HT=4
DEFAULT_MIN_HASHSIZE=1e6

def build_common_args():

    parser = argparse.ArgumentParser(description=
                                     'Build & load a counting Bloom filter.')

    env_ksize = os.environ.get('KHMER_KSIZE', DEFAULT_K)
    env_n_hashes = os.environ.get('KHMER_N_HASHES', DEFAULT_N_HT)
    env_hashsize = os.environ.get('KHMER_MIN_HASHSIZE', DEFAULT_MIN_HASHSIZE)

    parser.add_argument('-q', '--quiet', dest='quiet', default=False,
                        action='store_true')
    parser.add_argument('--ksize', '-k', type=int, dest='ksize',
                        default=env_ksize,
                        help='k-mer size to use')
    parser.add_argument('--n_hashes', '-N', type=int, dest='n_hashes',
                        default=env_n_hashes,
                        help='number of hash tables to use')
    parser.add_argument('--hashsize', '-x', type=float, dest='min_hashsize',
                        default=env_hashsize,
                        help='lower bound on hashsize to use')

    return parser

def parse_args(parser):
    args = parser.parse_args()

    if not args.quiet:
        if args.min_hashsize == DEFAULT_MIN_HASHSIZE:
            print>>sys.stderr, "** WARNING: hashsize is default!  You absodefly want to increase this!\n** Please read the docs!"

        print>>sys.stderr, '\nPARAMETERS:'
        print>>sys.stderr, ' - kmer size =    %d \t\t(-k)' % args.ksize
        print>>sys.stderr, ' - n hashes =     %d \t\t(-N)' % args.n_hashes
        print>>sys.stderr, ' - min hashsize = %-5.2g \t(-x)' % args.min_hashsize
        print>>sys.stderr, ''
        print>>sys.stderr, 'Estimated memory usage is %.2g bytes (n_hashes x min_hashsize)' % (args.n_hashes * args.min_hashsize)
        print>>sys.stderr, '-'*8

    return args

def main():
    parser = build_common_args()
    parser.add_argument('input_filename')
    parser.add_argument('output_filename')

    args = parse_args(parser)

    K=args.ksize
    HT_SIZE=args.min_hashsize
    N_HT=args.n_hashes

    input_name = args.input_filename
    output_name = args.output_filename

    print 'making hashtable'
    ht = khmer.new_counting_hash(K, HT_SIZE, N_HT)

    outfp = open(output_name, 'w')

    discarded = 0
    for n, record in enumerate(screed.open(input_name)):
        if n > 0 and n % 10000 == 0:
            print '...', n, discarded, int(discarded / float(n) * 100.)

        if len(record.sequence) < K:
            continue

        med, _, _ = ht.get_median_count(record.sequence)

        if med < DESIRED_COVERAGE:
            ht.consume(record.sequence)
            outfp.write('>%s\n%s\n' % (record.name, record.sequence))
        else:
            discarded += 1

    print 'consumed', n, 'discarded', discarded, 'percent', int(discarded / float(n) * 100.)


if __name__ == '__main__':
    main()
