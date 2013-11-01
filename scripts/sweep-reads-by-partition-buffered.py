#! /w/khmer_dev/bin/python

import screed
import sys
import argparse
import time
import khmer
from khmer.counting_args import build_construct_args, DEFAULT_MIN_HASHSIZE

# little class to store sequence information for the buffering class
class Seq:
    def __init__(self, name, color, seq):
        self.name = name
        self.color = color
        self.seq = seq
    def write(self, fp):
        fp.write('>{}\t{}\n{}\n'.format(self.name, self.color, self.seq))

# stores reads in memory and flushes them to their approriate files
# when certain criteria are met
# Basic idea is to buffer some number of reads in memory, then dump them all at once
# Hope that each file acrues, on average, BUFFER_SIZE / NUM_PARTS reads
# ie, if we buffer 1000000 reads, and we have 100000 partitions or colors,
# we should expect the mean buffer size to be 10 reads
class ReadBuffer:

    def __init__(self, max_buffers=10000, max_reads=1000000, est_files=100000, output_pref='reads_'):
        self.buffers = {}
        self.buffer_counts = {}
        self.max_buffers = max_buffers
        self.max_reads = max_reads

        self.est_files = est_files
        self.output_pref = output_pref
        self.buffer_flush = self.max_reads / self.est_files

        self.cur_reads = 0
        self.cur_files = 0

    def add_seq(self, seq):
        color = seq.color
        if color in self.buffers:
            count = self.buffer_counts[color]
            self.buffers[color].append(seq)
            self.buffer_counts[color] += 1
            if count > self.buffer_flush:
                self.flush_buffer(color)
                self.del_buffer(color)

        else:
            self.buffers[color] = [seq]
            self.buffer_counts[color] = 1
        self.cur_reads += 1
        if self.cur_reads > self.max_reads:
            self.flush_all()
        if len(self.buffers) > self.max_buffers:
            #self.clean_buffers(2)
            self.flush_all()
    
    def flush_buffer(self, color):
        with open('{}_{}.fa'.format(self.output_pref, color), 'a') as outfp:
            for read in self.buffers[color]:
                read.write(outfp)
                self.cur_reads -= 1
            
    def del_buffer(self, color):
        del self.buffer_counts[color]
        del self.buffers[color]

    def flush_all(self):
        print >>sys.stderr, '** reached max buffer size, flushing all to files...'
        for color in self.buffers:
            self.flush_buffer(color)
        colors = self.buffers.keys()
        for color in colors:
            self.del_buffer(color)
        del colors
        assert self.cur_reads == 0

    def clean_buffers(self, cutoff):
        print >>sys.stderr, '** flushing low-abundance buffers...'
        flushed = []
        for color in self.buffers:
            if self.buffer_counts[color] < cutoff:
                self.flush_buffer(color)
                flushed.append(color)
        for color in flushed:
            self.del_buffer(color)
        del flushed

def main():

    parser = build_construct_args()
    parser.add_argument('-i', '--input_fastp',dest='input_fastp')
    parser.add_argument('-r', '--traversal_range', type=int, dest='traversal_range')
    parser.add_argument('-b', '--buffer_size', dest='buffer_size', type=int)
    parser.add_argument('-e', '--files_estimate', dest='files_estimate', type=int)
    parser.add_argument('-o', '--output_prefix', dest='output_prefix')
    parser.add_argument('-m', '--max_buffers', dest='max_buffers', type=int)
    parser.add_argument('input_files', nargs='+')
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
    input_fastp = args.input_fastp

    max_buffers = args.max_buffers
    output_pref = args.output_prefix
    buf_size = args.buffer_size
    est = args.files_estimate
    input_files = args.input_files

    output_buffer = ReadBuffer(max_buffers=max_buffers, max_reads=buf_size, est_files=est, output_pref=output_pref)

	# file for multicolored reads, just keep this one around the whole time
    multi_fp = open('{}_multi.fp'.format(output_pref), 'a')
    orphaned_fp = open('{}_orphaned.fa'.format(output_pref), 'a')

	# consume the partitioned fasta with which to color the graph
    ht = khmer.new_hashbits(K, HT_SIZE, N_HT)
    print >>sys.stderr, 'consuming fastp...'
    ht.consume_partitioned_fasta_and_tag_with_colors(input_fastp)

    color_number_dist = []
    
    n_orphaned = 0
    n_colored = 0
    n_mcolored = 0
    try:
        total_t = time.clock()
        start_t = time.clock()
        for read_file in input_files:
            print >>sys.stderr,'** sweeping {read_file} for colors...'.format(read_file=read_file)
            file_t = 0.0
            for n, record in enumerate(screed.open(read_file)):

                if n % 50000 == 0:
                    end_t = time.clock()
                    batch_t = end_t - start_t
                    file_t += batch_t
                    print >>sys.stderr, '\tswept {n} reads [{nc} colored, {no} orphaned] ** {sec}s ({sect}s total)' \
                                        .format(n=n, nc=n_colored, no=n_orphaned, sec=batch_t, sect=file_t)
                    start_t = time.clock()
                seq = record.sequence
                name = record.name
                
                colors = ht.sweep_color_neighborhood(seq, traversal_range)
                color_number_dist.append(len(colors))
                if colors:
                    n_colored += 1
                    if len(colors) > 1:
                        multi_fp.write('>{}\t{}\n{}\n'.format(name, '\t'.join([str(c) for c in colors]), seq))
                    else:
                        output_buffer.add_seq(Seq(name, colors[0], seq))
                else:
                    n_orphaned += 1
                    orphaned_fp.write('>{}\n{}\n'.format(name, seq))

    except IOError as e:
        print >>sys.stderr, 'ERROR:', e
        print >>sys.stderr, '** exiting...'
    
	total_t = time.clock() - total_t

    print >>sys.stderr, 'swept {n_reads} for colors...'.format(n_reads=n)
    print >>sys.stderr, '...with {nc} colored and {no} orphaned'.format(
                                    nc=n_colored, no=n_orphaned)
    print >>sys.stderr, '...and {nmc} multicolored'.format(nmc=n_mcolored)
    
    print >>sys.stderr, '** outputting color number distribution...'
    with open('color_dist.txt', 'wb') as outfp:
        for nc in color_number_dist:
            outfp.write('{nc}\n'.format(nc=nc))
    
if __name__ == '__main__':
    main()
