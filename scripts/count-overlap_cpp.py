#! /usr/bin/env python
"""
Count the overlap k-mers, which are the k-mers apperaring in two sequence datase
ts.

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


###
DEFAULT_K=32
DEFAULT_N_HT=4
DEFAULT_HASHSIZE=1e6


def main():
    parser = argparse.ArgumentParser(description=
                                     'Use bloom filter to count overlap k-mers')
    env_ksize = os.environ.get('KHMER_KSIZE', DEFAULT_K)
    env_n_hashes = os.environ.get('KHMER_N_HASHES', DEFAULT_N_HT)
    env_hashsize = os.environ.get('KHMER_MIN_HASHSIZE', DEFAULT_HASHSIZE)
    parser.add_argument('-q', '--quiet', dest='quiet', default=False,
                        action='store_true')
    parser.add_argument('--ksize', '-k', type=int, dest='ksize',
                        default=env_ksize,
                        help='k-mer size to use')
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
            print>>sys.stderr, "** WARNING: hashsize is default!  You absodefly\
want to increase this!\n** Please read the docs!"
        print>>sys.stderr, '\nPARAMETERS:'
        print>>sys.stderr, ' - kmer size =    %d \t\t(-k)' % args.ksize
        print>>sys.stderr, ' - n hashes =     %d \t\t(-N)' % args.n_hashes
        print>>sys.stderr, ' - hashsize = %-5.2g \t(-x)' % args.hashsize
        print>>sys.stderr, 'Estimated memory usage is %.2g bytes (n_hashes x \
hashsize / 8)' % (args.n_hashes * args.hashsize / 8.)
        print>>sys.stderr, '-'*8

    K=args.ksize
    HT_SIZE=args.hashsize
    N_HT=args.n_hashes
    htfile = args.htfile
    fafile = args.fafile
    output_filename = args.report_filename

    print 'loading hashbits from', htfile
    ht2 = khmer.load_hashbits(htfile)
    K = ht2.ksize()

    output = open(output_filename, 'w')
    ht1 = khmer.new_hashbits(K, HT_SIZE, N_HT)
    tuple = ht1.count_overlap(fafile,ht2)
    print tuple
    output.write('unique k-mers in dataset2:'+str(tuple[0])+'\n')
    output.write('overlap k-mers:'+str(tuple[1]))





if __name__ == '__main__':
    main()
