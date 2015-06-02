#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
"""
@@
"""
from __future__ import print_function

import sys
import screed
import khmer
from khmer.khmer_args import build_counting_args, DEFAULT_MIN_TABLESIZE

DEFAULT_LOWER_CUTOFF = 2000
DEFAULT_UPPER_CUTOFF = 65535

###


def main():
    parser = build_construct_args()
    parser.add_argument('-l', '--lower-cutoff', type=int, dest='lower_cutoff',
                        default=DEFAULT_LOWER_CUTOFF)
    parser.add_argument('-u', '--upper-cutoff', type=int, dest='upper_cutoff',
                        default=DEFAULT_UPPER_CUTOFF)

    parser.add_argument('output_filename')
    parser.add_argument('input_filename')

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
        print('Estimated memory usage is %.2g bytes " \
            "(n_hashes x min_hashsize)' % (
            args.n_hashes * args.min_hashsize), file=sys.stderr)
        print('-' * 8, file=sys.stderr)

    K = args.ksize
    HT_SIZE = args.min_hashsize
    N_HT = args.n_hashes

    output = args.output_filename
    input = args.input_filename

    print('lower cutoff:', args.lower_cutoff)
    print('upper cutoff:', args.upper_cutoff)
    print('Saving stoptags to %s' % output)
    print('Loading sequences in %s' % input)

    ###

    print('making hashtable')
    ht = khmer.new_counting_hash(K, HT_SIZE, N_HT)
    ht.set_use_bigcount(True)

    print('consuming input', input)
    hb = ht.collect_high_abundance_kmers(input,
                                         args.lower_cutoff,
                                         args.upper_cutoff)

    print('saving stoptags', output)
    hb.save_stop_tags(output)

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
