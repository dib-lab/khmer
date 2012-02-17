#! /usr/bin/env python
"""
Count the overlap k-mers, which are the k-mers apperaring in two datasets.

usage: count-overlap.py [-h] [-q] [--ksize KSIZE] [--n_hashes N_HASHES]
                        [--hashsize HASHSIZE] [--curve CURVE_REPORT_FILENAME]
                        first_filename second_filename report_filename

Use '-h' for parameter help.

Note: "curve_report_filename" is optional. If it is provided, the script will
generate a report of the increase of overlap k-mers as the number of sequences
in the second dataset increases. The report can be used to generate a curve to
show that trend easily.

The report file contains the numbers of unique k-mers in the two datasets
seperately and the number of overlap k-mers appearing in both datasets..
"""
import khmer
import sys
import screed
from screed.fasta import fasta_iter
import argparse
import os
import math
    
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
    parser.add_argument('first_filename')
    parser.add_argument('second_filename')
    parser.add_argument('report_filename')
    parser.add_argument('--curve','-c',type=str,dest='curve_report_filename',
                        default = None,
                        help = 'report to generate the curve of the increase of'
                        +'overlap k-mers, optional, no such report to generate'
                        +'if no such argument is provided')
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
    filename = args.first_filename
    filename2 = args.second_filename
    file_result = args.report_filename
    file_curve = args.curve_report_filename
    count_overlap(K,HT_SIZE,N_HT,filename,filename2,file_result,file_curve)

def count_overlap(K, HT_SIZE, N_HT, filename, filename2, file_result,
                  file_curve):
    if file_curve:
        count = 0
        for n, record in enumerate(screed.open(filename2)):
            count = count + 1
        max_count = count / 100
        file3 = open(file_curve, 'w')

    # consume first dataset
    ht = khmer.new_hashbits(K, HT_SIZE, N_HT)
    n_unique = 0

    for n, record in enumerate(screed.open(filename)):
        sequence = record['sequence']
        seq_len = len(sequence)
        for n in range(0, seq_len+1-K):
            kmer = sequence[n:n+K]
            if (not ht.get(kmer)):
                n_unique += 1
            ht.count(kmer)

    print filename,'has been consumed.'
    fpr = (1 - math.exp(-n_unique/HT_SIZE))**N_HT
    printout1 = """\
%s:
# of unique k-mers: %d
# of occupied bin: %d
false positive rate: %e
""" % (filename, n_unique, ht.n_occupied(), fpr)

    # consume second dataset
    ht2 = khmer.new_hashbits(K, HT_SIZE, N_HT)
    n_unique = 0
    n_overlap = 0
    seq_count = 0

    for n, record in enumerate(screed.open(filename2)):
        sequence = record['sequence']
        seq_len = len(sequence)
        for n in range(0,seq_len+1-K):
            kmer = sequence[n:n+K]
            if (not ht2.get(kmer)):
                n_unique+=1
                if (ht.get(kmer)):
                    n_overlap+=1
            ht2.count(kmer)
        if file_curve:
            seq_count = seq_count + 1
            if seq_count == max_count:
                string = str(n_unique) + ' ' + str(n_overlap) + '\n'
                file3 = open(file_curve, 'a')
                file3.write(string)
                file3.close()
                seq_count = 0

    print filename2, 'has been consumed.'
    fpr = (1 - math.exp(-n_unique/HT_SIZE))**N_HT
    printout2 = """\
%s:
# of unique k-mers: %d
# of occupied bin: %d
false positive rate: %e
===============
# of overlap unique k-mers: %d""" % \
    (filename2, n_unique, ht2.n_occupied(), fpr, n_overlap)

    file_result_object = open(file_result,'w')
    file_result_object.write(printout1)
    file_result_object.write(printout2)

if __name__ == '__main__':
    main()


