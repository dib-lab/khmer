#! /usr/bin/env python
"""
get the number of reads with specific median k-mer abundance as reads are consumed

usage: get_median_curve.py [-h] [-q] [--ksize KSIZE] [--n_hashes N_HASHES]
                        [--hashsize HASHSIZE]   file_to_count report_file

Use '-h' for parameter help.

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
                                     'get the number of reads with specific median k-mer abundance as reads are consumed')
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
    parser.add_argument('filename')
    parser.add_argument('fileout')


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
    filename = args.filename
    fileout = args.fileout
    count_median(K,HT_SIZE,N_HT,filename,fileout)

def count_median(K,HT_SIZE,N_HT,filename,fileout):

    count = 0
    for n, record in enumerate(screed.open(filename)):
        count = count+1
    max_count = count/20
    print max_count

    ht = khmer.new_counting_hash(K, HT_SIZE, N_HT)
    ht.set_use_bigcount(True)
#    seq_array = []
    seq_count = 0
    median_array = [1,2,3,4,5]
    med={}
    
    for median in median_array:
        med[median] = 0
    #print med
    count = 0
    for n, record in enumerate(screed.open(filename)):
        sequence = record['sequence']
        ht.consume(sequence)
#        seq_array.append(sequence)
        seq_count = seq_count + 1
        if seq_count == max_count:
            count = count+1
            number_of_sequence_consumed = max_count*count
            counted_sequence = 0
            #print number_of_sequence_consumed
            for n2,record2 in enumerate(screed.open(filename)):
                counted_sequence = counted_sequence+1
                sequence2 = record2['sequence']
                #print sequence2
#for seq in seq_array:
                a, b, c = ht.get_median_count(sequence2)
                #print a,b,c
                for median in median_array:
                    if a == median:
                        #print "hit!"
                        med[a] = med[a]+1
                if counted_sequence == number_of_sequence_consumed:
                    break
                
            #print med
            fileout_obj = open(fileout,'a')
            print_line = str(number_of_sequence_consumed)
            for median in median_array:
                print_line = print_line+ '\t'+str(med[median])+'\t'
            print_line = print_line+'\n'
            fileout_obj.write(print_line)
            fileout_obj.close()
            seq_count = 0
            med={}
            for median in median_array:
                med[median] = 0

if __name__ == '__main__':
    main()


