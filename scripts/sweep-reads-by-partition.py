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
from khmer.counting_args import build_construct_args, DEFAULT_MIN_HASHSIZE

DEFAULT_PPF = 1

def write_read(fp, seq, name, color):
    fp.write('>{name}\t{color}\n{seq}\n'.format(seq=seq, name=name, color=color))

def main():
    parser = build_construct_args()
    parser.add_argument('-p', '--partitions_per_file', 
                        dest='partitions_per_file', default=DEFAULT_PPF)
    parser.add_argument('-i', '--input_fastp', dest='input_fastp')
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
            (n_hashes x min_hashsize)'.format(prod=args.n_hashes*args.min_hashsize)
        print >>sys.stderr, '-' * 8
    
    K = args.ksize
    HT_SIZE = args.min_hashsize
    N_HT = args.n_hashes
    
    input_reads = args.input_reads
    input_fastp = args.input_fastp
    ppf = args.partitions_per_file
    
    ht = khmer.new_hashbits(K, HT_SIZE, N_HT)
    ht.consume_partitioned_fasta_and_tag_with_colors(input_fastp)
    
    cur_colors = []
    color_to_fp_dict = {}
    cur_fp = file
    
    color_number_dist = []
    
    n_orphaned = 0
    n_colored = 0
    n_mcolored = 0
    n_files = 0
    try:
        for read_file in input_reads:
            print >>sys.stderr,'** sweeping {read_file} for colors...'.format(read_file=read_file)
            
            for n, record in enumerate(screed.open(read_file)):
                if n % 50000 == 0:
                    print >>sys.stderr, '\tswept {n} reads [{nc} colored, {no} orphaned]' \
                                        .format(n=n, nc=n_colored, no=n_orphaned)
                seq = record.sequence
                name = record.name
                
                colors = ht.sweep_color_neighborhood(seq)
                color_number_dist.append(len(colors))
                if colors:
                    n_colored += 1
                    if len(colors) > 1:
                        n_mcolored += 1
                    for color in colors:
                        # do we have a file for this color already? use it!
                        if color in color_to_fp_dict:
                            fp = color_to_fp_dict[color]
                            write_read(fp, seq, name, color)
                        # no file yet? make a new one
                        else:
                            if len(cur_colors) == 0:
                                #print '** opening new file...'
                                cur_fp = open('colored_reads_{fn}.fa'.format(fn=n_files),
                                              'wb')
                                              
                            color_to_fp_dict[color] = cur_fp
                            cur_colors.append(color)
                            write_read(cur_fp, seq, name, color)
                            n_files += 1
                            
                            if len(cur_colors) == ppf:
                                cur_colors = []
                else:
                    n_orphaned += 1
            
        for key in color_to_fp_dict:
            if color_to_fp_dict[key]:
                color_to_fp_dict[key].close()

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
