#! /usr/bin/python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#

"""
Find all reads connected to the given contigs on a per-partition basis.

% python scripts/normalize-by-median.py -r <range> -i <contigs fastp> \
<reads1> <reads2> ... <readsN>

This script is very lenient on IO errors, due to the large number of file
operations needed. Thus, errors opening a file for buffer flush or writeing
a read to a file will not crash the program; instead, if there were errors,
the user will be warned at the end of execution. Errors with opening read files
are also handled -- we move on to the next read file if there is an error opening.

"""

import screed
import sys
import os
import argparse
import time
import khmer
from khmer.counting_args import build_construct_args, DEFAULT_MIN_HASHSIZE

DEFAULT_NUM_BUFFERS=50000
DEFAULT_BUFFER_SIZE=1000000
DEFAULT_NUM_PARTITIONS=100000
DEFAULT_OUT_PREF='reads_'
DEFAULT_RANGE=-1

MIN_HSIZE=4e7
MIN_KSIZE=21

# little class to store sequence information for the buffering class
class Seq:
    def __init__(self, name, color, seq):
        self.name = name
        self.color = color
        self.seq = seq
    
    def __repr__(self):
        return '''>{name}\t{color}\n
{seq}\n'''.format(name=self.name, color=self.color, seq=self.seq)

    def write(self, fp):
        try:
            fp.write('>{}\t{}\n{}\n'.format(self.name, self.color, self.seq))
        except IOError:
            print >>sys.stderr, 'Error writing {seq} to {fn}'.format(seq=self, fn=fp)
            return 1
        else:
            return 0

# stores reads in memory and flushes them to their approriate files
# when certain criteria are met
# Basic idea is to buffer some number of reads in memory, then dump them all at once
# Hope that each file acrues, on average, BUFFER_SIZE / NUM_PARTS reads
# ie, if we buffer 1000000 reads, and we have 100000 partitions or colors,
# we should expect the mean buffer size to be 10 reads
class ReadBuffer:

    def __init__(self, max_buffers, max_size, est_files, output_pref, outdir):
        self.buffers = {}
        self.buffer_counts = {}
        self.max_buffers = max_buffers
        self.max_size = max_size

        self.est_files = est_files
        self.output_pref = output_pref
        self.outdir = outdir
        self.buffer_flush = self.max_size / self.est_files

        self.cur_reads = 0
        self.cur_files = 0

        self.num_write_errors = 0
        self.num_file_errors = 0

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
        if self.cur_reads > self.max_size:
            self.flush_all()
        if len(self.buffers) > self.max_buffers:
            #self.clean_buffers(2)
            self.flush_all()
    
    def flush_buffer(self, color):
        fn = '{}_{}.fa'.format(self.output_pref, color)
        fpath = os.path.join(self.outdir, fn)
        try:
            outfp = open(fpath, 'a')
        except IOError as e:
            print >>sys.stderr, 'ERROR: {e}'.format(e=e)
            print >>sys.stderr, '*** Failed to open {fn} for buffer flush'.format(fpath)
            self.num_file_errors += 1
        else:
            for read in self.buffers[color]:
                self.num_write_errors += read.write(outfp)
                self.cur_reads -= 1
            outfp.close()

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
    parser.add_argument('-r', '--traversal_range', type=int, dest='traversal_range', \
                        default=DEFAULT_RANGE)
    parser.add_argument('-b', '--buffer_size', dest='buffer_size', type=int, \
                        default=DEFAULT_BUFFER_SIZE)
    parser.add_argument('-e', '--files_estimate', dest='files_estimate', type=int, \
                        default=DEFAULT_NUM_PARTITIONS)
    parser.add_argument('-o', '--output_prefix', dest='output_prefix',
                        default=DEFAULT_OUT_PREF)
    parser.add_argument('-m', '--max_buffers', dest='max_buffers', type=int, \
                        default=DEFAULT_NUM_BUFFERS)
    parser.add_argument('-d', '--debug', dest='debug', default=None)
    parser.add_argument('input_files', nargs='+')
    args = parser.parse_args()
    
   
    K = args.ksize
    HT_SIZE = args.min_hashsize
    if HT_SIZE < MIN_HSIZE:
        HT_SIZE = MIN_HSIZE
    if K < MIN_KSIZE:
        K = MIN_KSIZE
    N_HT = args.n_hashes

    if not args.quiet:
        if args.min_hashsize == DEFAULT_MIN_HASHSIZE:
            print >>sys.stderr, \
                "** WARNING: hashsize is default!  " \
                "You absodefly want to increase this!\n** " \
                "Please read the docs!"

        print >>sys.stderr, '\nPARAMETERS:'
        print >>sys.stderr, \
            ' - kmer size =    {ksize:d} \t\t(-k)'.format(ksize=K)
        print >>sys.stderr, \
            ' - n hashes =     {nhash:d} \t\t(-N)'.format(nhash=args.n_hashes)
        print >>sys.stderr, \
            ' - min hashsize = {mh:-5.2g} \t(-x)'.format(mh=HT_SIZE)
        print >>sys.stderr, ''
        print >>sys.stderr, \
            'Estimated memory usage is {prod:.2g} bytes \
            (n_hashes x min_hashsize / 8)'.format(prod=args.n_hashes*args.min_hashsize/8)
        print >>sys.stderr, '-' * 8
    
    traversal_range = args.traversal_range
    input_fastp = args.input_fastp
    outdir = os.path.dirname(input_fastp)

    max_buffers = args.max_buffers
    output_pref = args.output_prefix
    buf_size = args.buffer_size
    est = args.files_estimate
    input_files = args.input_files

    debug = args.debug
    if debug:
        import yep

    output_buffer = ReadBuffer(max_buffers, buf_size, est, output_pref, outdir)

	# file for multicolored reads, just keep this one around the whole time
    multi_fn = os.path.join(outdir, '{}_multi.fa'.format(output_pref))
    try:
        multi_fp = open(multi_fn, 'a')
    except IOError as e:
        print >>sys.stderr, 'ERROR: {e}'.format(e=e)
        print >>sys.stderr, '*** Failed to open {fn}'.format(multi_fn)
    orphaned_fn = os.path.join(outdir, '{}_orphaned.fa'.format(output_pref))
    try:
        orphaned_fp = open(orphaned_fn, 'a')
    except IOError as e:
        print >>sys.stderr, 'ERROR: {e}'.format(e=e)
        print >>sys.stderr, '*** Failed to open {fn}'.format(orphaned_fn)

	# consume the partitioned fasta with which to color the graph
    ht = khmer.new_hashbits(K, HT_SIZE, N_HT)
    print >>sys.stderr, 'consuming fastp...'
    if debug:
        yep.start(debug)
    ht.consume_partitioned_fasta_and_tag_with_colors(input_fastp)

    color_number_dist = []
    
    n_orphaned = 0
    n_colored = 0
    n_mcolored = 0

    total_t = time.clock()
    start_t = time.clock()
    for read_file in input_files:
        print >>sys.stderr,'** sweeping {read_file} for colors...'.format(read_file=read_file)
        file_t = 0.0
        try:
            read_fp = screed.open(read_file)
        except IOError as e:
            print >>sys.stderr, 'ERROR:', e
            print >>sys.stderr, '*** Could not open {fn}, skipping...'.format(fn=read_file)
        else:
            for n, record in enumerate(read_fp):
                if n % 50000 == 0:
                    end_t = time.clock()
                    batch_t = end_t - start_t
                    file_t += batch_t
                    print >>sys.stderr, '\tswept {n} reads [{nc} colored, {no} orphaned] \
                                        ** {sec}s ({sect}s total)' \
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
            output_buffer.flush_all()
            read_fp.close()

    # gotta output anything left in the buffers at the end!
    output_buffer.flush_all() 
    total_t = time.clock() - total_t

    multi_fp.close()
    orphaned_fp.close()
    if debug:
        yep.stop()
    if output_buffer.num_write_errors > 0 or output_buffer.num_file_errors > 0:
        print >>sys.stderr, 'WARNING: Sweep finished with errors!'
        print >>sys.stderr, '** {writee} reads not written'.format(writee=output_buffer.num_write_errors)
        print >>sys.stderr, '** {filee} errors opening files'.format(filee=output_buffer.num_file_errors)

    print >>sys.stderr, 'swept {n_reads} for colors...'.format(n_reads=n_colored+n_mcolored+n_orphaned)
    print >>sys.stderr, '...with {nc} colored and {no} orphaned'.format(
                                    nc=n_colored, no=n_orphaned)
    print >>sys.stderr, '...and {nmc} multicolored'.format(nmc=n_mcolored)
    
    print >>sys.stderr, '** outputting color number distribution...'
    with open('color_dist.txt', 'wb') as outfp:
        for nc in color_number_dist:
            outfp.write('{nc}\n'.format(nc=nc))
    
if __name__ == '__main__':
    main()
