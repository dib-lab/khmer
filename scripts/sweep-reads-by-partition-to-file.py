#! /w/khmer_dev/bin/python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
"""
Tag and color the given partitioned fasta, then find all reads in the neighborhood
of each partition and output to a file

% python scripts/normalize-by-median.py [ -p <partitions/file> ] -i <fastp> <reads1> <reads2> ...

Use '-h' for parameter help.
"""

import khmer
import screed
import sys
import time
from khmer.counting_args import build_construct_args, DEFAULT_MIN_HASHSIZE

MAX_FILES=512
READS_PER_FILE = 100000000

def write_read(fp, seq, name, color):
    fp.write('>{name}\t{color}\n{seq}\n'.format(seq=seq, name=name, color=color))

def main():
    parser = build_construct_args()
    #parser.add_argument('-p', '--partitions_per_file', 
    #                    dest='partitions_per_file', default=DEFAULT_PPF)
    parser.add_argument('-i', '--input_fastp',dest='input_fastp')
    parser.add_argument('-r', '--traversal_range', type=int, dest='traversal_range')
    parser.add_argument('input_reads', nargs='+')
    args = parser.parse_args()
    
    if not args.quiet:
        if args.min_hashsize == DEFAULT_MIN_HASHSIZE:
            print >>sys.stderr, \
                "** WARNING: hashsize is default!  " \
                "You absodefly want to increase this!\n** " \
                "Please read the docs!"

        print >>sys.stderr, '\nPARAMETERS:'
        print >>sys.stderr, \
            ' - kmer size =    {ksize:d} \t\t(-k)'.format(ksize=args.ksize)
        print >>sys.stderr, \
            ' - n hashes =     {nhash:d} \t\t(-N)'.format(nhash=args.n_hashes)
        print >>sys.stderr, \
            ' - min hashsize = {mh:-5.2g} \t(-x)'.format(mh=args.min_hashsize)
        print >>sys.stderr, ''
        print >>sys.stderr, \
            'Estimated memory usage is {prod:.2g} bytes \
            (n_hashes x min_hashsize / 8)'.format(prod=args.n_hashes*args.min_hashsize/8)
        print >>sys.stderr, '-' * 8
    
    K = args.ksize
    HT_SIZE = args.min_hashsize
    N_HT = args.n_hashes
    
    traversal_range = args.traversal_range
    input_reads = args.input_reads
    input_fastp = args.input_fastp
    
    ht = khmer.new_hashbits(K, HT_SIZE, N_HT)
    print >>sys.stderr, 'consuming fastp...'
    ht.consume_partitioned_fasta_and_tag_with_colors(input_fastp)
 
    color_number_dist = []
    
    n_orphaned = 0
    n_colored = 0
    n_mcolored = 0
    n_files = 0
    try:
        outfp = open('colored_reads_0.fa', 'wb')
        start_t = time.clock()
        for read_file in input_reads:
            print >>sys.stderr,'** sweeping {read_file} for colors...'.format(read_file=read_file)
            total_t = 0.0
            for n, record in enumerate(screed.open(read_file)):
                if n % 50000 == 0:
                    end_t = time.clock()
                    batch_t = end_t - start_t
                    total_t += batch_t
                    print >>sys.stderr, '\tswept {n} reads [{nc} colored, {no} orphaned] ** {sec}s ({sect}s total)' \
                                        .format(n=n, nc=n_colored, no=n_orphaned, sec=batch_t, sect=total_t)
                    start_t = time.clock()
                seq = record.sequence
                name = record.name
                
                colors = ht.sweep_color_neighborhood(seq, traversal_range)
                color_number_dist.append(len(colors))
                if colors:
                    n_colored += 1
                    if len(colors) > 1:
                        n_mcolored += 1
                    for color in colors:
                        write_read(outfp, seq, name, color)
                else:
                    n_orphaned += 1

                if n_colored % READS_PER_FILE == 0 and n_colored != 0:
                    n_files += 1
                    outfp = open('colored_reads_{}.fa'.format(n_files), 'wb')

    except IOError as e:
        print >>sys.stderr, 'ERROR:', e
        print >>sys.stderr, '** exiting...'
        
    print >>sys.stderr, 'swept {n_reads} for colors...'.format(n_reads=n)
    print >>sys.stderr, '...with {nc} colored and {no} orphaned'.format(
                                    nc=n_colored, no=n_orphaned)
    print >>sys.stderr, '...and {nmc} multicolored'.format(nmc=n_mcolored)
    print >>sys.stderr, '...to {nf} files'.format(nf=n_files)
    
    print >>sys.stderr, '** outputting color number distribution...'
    with open('color_dist.txt', 'wb') as outfp:
        for nc in color_number_dist:
            outfp.write('{nc}\n'.format(nc=nc))
    
if __name__ == '__main__':
    main()
