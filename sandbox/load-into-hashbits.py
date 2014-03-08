#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
"""
Build a graph Bloom filter from the given sequences, save in <htname>.

% python scripts/load-into-hashbits.py <htname> <data1> [ <data2> <...> ]

Parameters to adjust: K, HT_SIZE.  HT_SIZE should be set to about 2x the
available system memory.
"""

import sys
import screed
import os
import khmer
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
        if args.min_hashsize == DEFAULT_MIN_HASHSIZE:
            print >>sys.stderr, "** WARNING: hashsize is default!  " \
                "You absodefly want to increase this!\n** " \
                "Please read the docs!"

        print >>sys.stderr, '\nPARAMETERS:'
        print >>sys.stderr, ' - kmer size =    %d \t\t(-k)' % args.ksize
        print >>sys.stderr, ' - n hashes =     %d \t\t(-N)' % args.n_hashes
        print >>sys.stderr, ' - min hashsize = %-5.2g \t(-x)' % \
            args.min_hashsize
        print >>sys.stderr, ''
        print >>sys.stderr, 'Estimated memory usage is %.2g bytes ' \
            '(n_hashes x min_hashsize / 8 bits/byte)' % (
            args.n_hashes * args.min_hashsize / 8.)
        print >>sys.stderr, '-' * 8

    return args

###


def main():
    parser = build_common_args()
    parser.add_argument('output_filename')
    parser.add_argument('input_filenames', nargs='+')

    args = parse_args(parser)

    K = args.ksize
    HT_SIZE = args.min_hashsize
    N_HT = args.n_hashes

    base = args.output_filename
    filenames = args.input_filenames

    print 'Saving hashtable to %s' % base
    print 'Loading kmers from sequences in %s' % repr(filenames)

    ###

    print 'making hashtable'
    ht = khmer.new_hashbits(K, HT_SIZE, N_HT)

    for n, filename in enumerate(filenames):
        print 'consuming input', filename
        ht.consume_fasta(filename)

        if n > 0 and n % 10 == 0:
            print 'mid-save', base
            ht.save(base)
            open(base + '.info', 'w').write('through %s' % filename)

    print 'saving', base
    ht.save(base)
    open(base + '.info', 'w').write('through end: %s' % filename)

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
