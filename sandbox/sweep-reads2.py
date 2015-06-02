#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
"""
Use a set of query reads to sweep out overlapping reads from another file.

% python scripts/sweep-reads2.py <query reads> <search reads>

Results end up in <search reads>.sweep2.

Use '-h' for parameter help.
"""
from __future__ import print_function

import sys
import khmer
import os.path
import screed
from khmer.khmer_args import (build_hashbits_args, DEFAULT_MIN_TABLESIZE)


def main():
    parser = build_construct_args()
    parser.add_argument('input_filename')
    parser.add_argument('read_filename')

    args = parser.parse_args()

    if not args.quiet:
        if args.min_hashsize == DEFAULT_MIN_HASHSIZE:
            print("** WARNING: hashsize is default!  " \
                "You absodefly want to increase this!\n** " \
                "Please read the docs!", file=sys.stderr)

        print('\nPARAMETERS:', file=sys.stderr)
        print(' - kmer size =    %d \t\t(-k)' % args.ksize, file=sys.stderr)
        print(' - n hashes =     %d \t\t(-N)' % args.n_hashes, file=sys.stderr)
        print(' - min hashsize = %-5.2g \t(-x)' % \
            args.min_hashsize, file=sys.stderr)
        print('', file=sys.stderr)
        print('Estimated memory usage is %.2g bytes ' \
            '(n_hashes x min_hashsize / 8)' % (
                args.n_hashes * args.min_hashsize / 8.), file=sys.stderr)
        print('-' * 8, file=sys.stderr)

    K = args.ksize
    HT_SIZE = args.min_hashsize
    N_HT = args.n_hashes

    inp = args.input_filename
    readsfile = args.read_filename

    outfile = os.path.basename(readsfile) + '.sweep2'
    outfp = open(outfile, 'w')

    # create a hashbits data structure
    ht = khmer.new_hashbits(K, HT_SIZE, N_HT)

    # load contigs, connect into N partitions
    print('loading input reads from', inp)
    ht.consume_fasta(inp)

    print('starting sweep.')

    n = 0
    m = 0
    for record in screed.open(readsfile):
        if len(record.sequence) < K:
            continue

        if n % 10000 == 0:
            print('...', n, m)

        count = ht.get_median_count(record.sequence)[0]
        if count:
            m += 1
            outfp.write('>%s\n%s\n' % (record.name, record.sequence))
        n += 1

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
