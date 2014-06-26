#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
"""
Eliminate reads with median k-mer abundance higher than DESIRED_COVERAGE.

Parameters to adjust: DESIRED_COVERAGE,
"""

import sys
import screed
import os
import khmer
import random

DEFAULT_DESIRED_COVERAGE = 20

import argparse

DEFAULT_K = 32
DEFAULT_N_HT = 4
DEFAULT_MIN_HASHSIZE = 1e6


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
        if args.min_hashsize == DEFAULT_MIN_HASHSIZE and not args.loadhash:
            print >>sys.stderr, \
                "** WARNING: hashsize is default!  ", \
                "You absodefly want to increase this!\n** " \
                "Please read the docs!"

        print >>sys.stderr, '\nPARAMETERS:'
        print >>sys.stderr, ' - kmer size =    %d \t\t(-k)' % args.ksize
        print >>sys.stderr, ' - n hashes =     %d \t\t(-N)' % args.n_hashes
        print >>sys.stderr, ' - min hashsize = %-5.2g \t(-x)' % \
            args.min_hashsize
        print >>sys.stderr, ''
        print >>sys.stderr, \
            'Estimated memory usage is %.2g bytes ' \
            '(n_hashes x min_hashsize)' % (
                args.n_hashes * args.min_hashsize)
        print >>sys.stderr, '-' * 8

    return args


def main():
    parser = build_common_args()
    parser.add_argument('input_filenames', nargs='+')
    parser.add_argument('-C', '--cutoff', type=int, dest='cutoff',
                        default=DEFAULT_DESIRED_COVERAGE)
    parser.add_argument('-s', '--savehash', dest='savehash', default='')
    parser.add_argument('-l', '--loadhash', dest='loadhash',
                        default='')

    args = parse_args(parser)

    K = args.ksize
    HT_SIZE = args.min_hashsize
    N_HT = args.n_hashes
    DESIRED_COVERAGE = args.cutoff

    input_name_list = args.input_filenames

    if args.loadhash:
        print 'loading hashtable from', args.loadhash
        ht = khmer.load_counting_hash(args.loadhash)
    else:
        print 'making hashtable'
        ht = khmer.new_counting_hash(K, HT_SIZE, N_HT)

    total = 0
    discarded = 0
    for input_filename in input_name_list:
        output_name = os.path.basename(input_filename) + '.keep'
        outfp = open(output_name, 'w')

        for n, record in enumerate(screed.open(input_filename)):
            if n > 0 and n % 10000 == 0:
                print '... kept', total - discarded, 'of', total, ', or', \
                    int(100. - discarded / float(total) * 100.), '%'
                print '... in file', input_filename

            total += 1

            if len(record.sequence) < K:
                continue

            seq = record.sequence.replace('N', 'A')

            med, _, _ = ht.get_median_count(seq)

            if med < DESIRED_COVERAGE:
                ht.consume(seq)
                outfp.write('>%s\n%s\n' % (record.name, record.sequence))
            else:
                discarded += 1

        print 'DONE with', input_filename, '; kept', total - discarded, 'of', \
            total, 'or', int(100. - discarded / float(total) * 100.), '%'

    if args.savehash:
        print 'Saving hashfile through', input_filename
        print '...saving to', args.savehash
        ht.save(os.path.basename(args.savehash))

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
