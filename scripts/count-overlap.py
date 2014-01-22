#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
"""
Count the overlap k-mers, which are the k-mers apperaring in two sequence
datasets.

usage: count-overlap_cpp.py [-h] [-q] [--ksize KSIZE] [--n_hashes N_HASHES]
                        [--hashsize HASHSIZE]
                           1st_dataset(htfile) 2nd_dataset(fastafile) result

Use '-h' for parameter help.

"""
import khmer
import sys
import screed
from screed.fasta import fasta_iter
import argparse
import os
import math


#
DEFAULT_K = 32
DEFAULT_N_HT = 4
DEFAULT_HASHSIZE = 1e6


def main():
    parser = argparse.ArgumentParser(
        description='Use bloom filter to count overlap k-mers',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    env_ksize = os.environ.get('KHMER_KSIZE', DEFAULT_K)
    env_n_hashes = os.environ.get('KHMER_N_HASHES', DEFAULT_N_HT)
    env_hashsize = os.environ.get('KHMER_MIN_HASHSIZE', DEFAULT_HASHSIZE)
    parser.add_argument('-q', '--quiet', dest='quiet', default=False,
                        action='store_true')
    parser.add_argument('--ksize', '-k', type=int, dest='ksize',
                        default=env_ksize,
                        help='k-mer size to use '
                             '(should be the same as in htfile)')
    parser.add_argument('--n_hashes', '-N', type=int, dest='n_hashes',
                        default=env_n_hashes,
                        help='number of hash tables to use')
    parser.add_argument('--hashsize', '-x', type=float, dest='hashsize',
                        default=env_hashsize,
                        help='hashsize to use')
    parser.add_argument('htfile')
    parser.add_argument('fafile')
    parser.add_argument('report_filename')
    args = parser.parse_args()
    if not args.quiet:
        if args.hashsize == DEFAULT_HASHSIZE:
            print >>sys.stderr, \
                "** WARNING: hashsize is default!  " \
                "You absodefly want to increase this!\n** " \
                "Please read the docs!"
        print >>sys.stderr, '\nPARAMETERS:'
        print >>sys.stderr, ' - kmer size =    %d \t\t(-k)' % args.ksize
        print >>sys.stderr, ' - n hashes =     %d \t\t(-N)' % args.n_hashes
        print >>sys.stderr, ' - hashsize = %-5.2g \t(-x)' % args.hashsize
        print >>sys.stderr, \
            'Estimated memory usage is %.2g bytes (n_hashes x hashsize / 8)' \
            % (args.n_hashes * args.hashsize / 8.)
        print >>sys.stderr, '-' * 8

    K = args.ksize
    HT_SIZE = args.hashsize
    N_HT = args.n_hashes
    htfile = args.htfile
    fafile = args.fafile
    output_filename = args.report_filename
    curve_filename = output_filename + '.curve'

    print 'loading hashbits from', htfile
    ht1 = khmer.load_hashbits(htfile)
    K = ht1.ksize()

    output = open(output_filename, 'w')
    f_curve_obj = open(curve_filename, 'w')

    ht2 = khmer.new_hashbits(K, HT_SIZE, N_HT)

    (n_unique, n_overlap, list) = ht2.count_overlap(fafile, ht1)

    printout1 = """\
dataset1(ht file): %s
dataset2: %s

# of unique k-mers in dataset2: %d
# of overlap unique k-mers: %d

""" % (htfile, fafile, n_unique, n_overlap)
    output.write(printout1)

    figure_list1 = []
    figure_list2 = []

    for i in range(100):
        to_print = str(list[100 + i]) + ' ' + str(list[i]) + '\n'
        f_curve_obj.write(to_print)


if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
