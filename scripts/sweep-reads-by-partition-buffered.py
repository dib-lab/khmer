#!/usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#

"""
Find all reads connected to the given contigs on a per-partition basis.

% sweep-reads-by-partition.py -r <range> <contigs fastp> \
<reads1> <reads2> ... <readsN>

This script is very lenient on IO errors, due to the large number of file
operations needed. Thus, errors opening a file for buffer flush or writing
a read to a file will not crash the program; instead, if there were errors,
the user will be warned at the end of execution. Errors with opening read files
are also handled -- we move on to the next read file if there is an error
opening.
"""

import screed
import sys
import os
import time
import khmer
from khmer.counting_args import build_construct_args, DEFAULT_MIN_HASHSIZE


DEFAULT_NUM_BUFFERS = 50000
DEFAULT_MAX_READS = 1000000
DEFAULT_BUFFER_SIZE = 10
DEFAULT_OUT_PREF = 'reads_'
DEFAULT_RANGE = -1

MIN_HSIZE = 4e7
MIN_KSIZE = 21


def fmt_fasta(name, seq, labels=[]):
    return '>{name}\t{labels}\n{seq}\n'.format(name=name,
            labels='\t'.join([str(l) for l in labels]), seq=seq)


def write_seq(fp, name, seq, labels=[]):
    try:
        fp.write(fmt_fasta(name, seq, labels=labels))
    except IOError:
        print >>sys.stderr, 'Error writing {read}'.format(
            read=fmt_fasta(name, seq, labels=labels))
        return 1
    else:
        return 0


class ReadBuffer:

    def __init__(self):
        self.buf = []

    def push(self, seq_str):
        self.buf.append(seq_str)

    def flush(self):
        return ''.join(self.buf)

    def is_full(self, full):
        if len(self.buf) >= full:
            return True
        else:
            return False

    def __len__(self):
        return len(self.buf)


class ReadBufferManager:

    def __init__(self, max_buffers, max_reads, max_size, output_pref, outdir):
        self.buffers = {}
        self.buffer_counts = {}
        self.max_buffers = max_buffers
        self.max_reads = max_reads

        self.output_pref = output_pref
        self.outdir = outdir
        self.buffer_flush = max_size

        self.cur_reads = 0
        self.cur_files = 0

        self.num_write_errors = 0
        self.num_file_errors = 0

        print >>sys.stderr, '''Init new ReadBuffer [
        Max Buffers: {num_bufs}
        Max Reads: {max_reads}
        Buffer flush: {buf_flush}
        ]'''.format(num_bufs=self.max_buffers, max_reads=self.max_reads,
                    buf_flush=self.buffer_flush)

    def flush_buffer(self, buf_id):
        fn = '{}_{}.fa'.format(self.output_pref, buf_id)
        fpath = os.path.join(self.outdir, fn)
        try:
            outfp = open(fpath, 'a')
        except IOError as e:
            print >>sys.stderr, '!! ERROR: {e} !!'.format(e=e)
            print >>sys.stderr, '*** Failed to open {fn} for buffer flush'.format(fn=fpath)
            self.num_file_errors += 1
        else:
            buf = self.buffers[buf_id]
            outfp.write(buf.flush())
            self.cur_reads -= len(buf)
            outfp.close()
            del self.buffers[buf_id]

    def queue(self, seq_str, buf_id):
        if buf_id in self.buffers:
            self.buffers[buf_id].push(seq_str)
            if self.buffers[buf_id].is_full(self.buffer_flush):
                self.flush_buffer(buf_id)
        else:
            new_buf = ReadBuffer()
            new_buf.push(seq_str)
            self.buffers[buf_id] = new_buf

        self.cur_reads += 1
        if self.cur_reads > self.max_reads:
            print >>sys.stderr, '** Reached max num reads...'
            self.flush_all()
        if len(self.buffers) > self.max_buffers:
            # self.clean_buffers(2)
            print >>sys.stderr, '** Reached max num buffers...'
            self.flush_all()

    def flush_all(self):
        print >>sys.stderr, '*** Flushing all to files...'
        buf_ids = self.buffers.keys()
        for buf_id in buf_ids:
            self.flush_buffer(buf_id)
        assert self.cur_reads == 0


def main():

    parser = build_construct_args('Takes a partitioned reference file \
                                  and a list of reads, and sorts reads \
                                  by which partition they connect to')
    parser.add_argument(
        '-r', '--traversal_range', type=int, dest='traversal_range',
        default=DEFAULT_RANGE)
    parser.add_argument('-b', '--buffer_size', dest='max_reads', type=int,
                        default=DEFAULT_MAX_READS)
    parser.add_argument('-l', '--buffer_length', dest='buffer_size', type=int,
                        default=DEFAULT_BUFFER_SIZE)
    parser.add_argument('-o', '--output_prefix', dest='output_prefix',
                        default=DEFAULT_OUT_PREF)
    parser.add_argument('-m', '--max_buffers', dest='max_buffers', type=int,
                        default=DEFAULT_NUM_BUFFERS)
    parser.add_argument(dest='input_fastp')
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
            (n_hashes x min_hashsize / 8)'.format(prod=args.n_hashes * HT_SIZE / 8)
        print >>sys.stderr, '-' * 8

    traversal_range = args.traversal_range
    input_fastp = args.input_fastp
    outdir = os.path.dirname(input_fastp)

    max_buffers = args.max_buffers
    output_pref = args.output_prefix
    buf_size = args.buffer_size
    max_reads = args.max_reads

    input_files = args.input_files

    output_buffer = ReadBufferManager(
        max_buffers, max_reads, buf_size, output_pref, outdir)

    # consume the partitioned fasta with which to label the graph
    ht = khmer.LabelHash(K, HT_SIZE, N_HT)
    print >>sys.stderr, 'consuming fastp...'
    ht.consume_partitioned_fasta_and_tag_with_labels(input_fastp)

    label_number_dist = []

    n_orphaned = 0
    n_labeled = 0
    n_mlabeled = 0

    total_t = time.clock()
    start_t = time.clock()
    for read_file in input_files:
        print >>sys.stderr, '** sweeping {read_file} for labels...'.format(
                                                        read_file=read_file)
        file_t = 0.0
        try:
            read_fp = screed.open(read_file)
        except IOError as e:
            print >>sys.stderr, '!! ERROR: !!', e
            print >>sys.stderr, '*** Could not open {fn}, skipping...'.format(
                                                                 fn=read_file)
        else:
            for n, record in enumerate(read_fp):
                if n % 50000 == 0:
                    end_t = time.clock()
                    batch_t = end_t - start_t
                    file_t += batch_t
                    print >>sys.stderr, '\tswept {n} reads [{nc} labeled, \
                                         {no} orphaned] \
                                        ** {sec}s ({sect}s total)' \
                                        .format(n=n, nc=n_labeled, 
                                                no=n_orphaned, 
                                                sec=batch_t, sect=file_t)
                    start_t = time.clock()
                seq = record.sequence
                name = record.name
                try:
                    labels = ht.sweep_label_neighborhood(seq, traversal_range)
                except ValueError as e:
                    print >>sys.stderr, '!! ERROR: {e} !!'.format(e=e)
                    print >>sys.stderr, 'Read length less than k-mer size'
                else:
                    seq_str = fmt_fasta(name, seq, labels)
                    label_number_dist.append(len(labels))
                    if labels:
                        n_labeled += 1
                        if len(labels) > 1:
                            output_buffer.queue(seq_str, 'multi')
                            n_mlabeled += 1
                        else:
                            output_buffer.queue(seq_str, labels[0])
                    else:
                        n_orphaned += 1
                        output_buffer.queue(seq_str, 'orphaned')
            print >>sys.stderr, '** End of file {fn}...'.format(fn=read_file)
            output_buffer.flush_all()
            read_fp.close()

    # gotta output anything left in the buffers at the end!
    print >>sys.stderr, '** End of run...'
    output_buffer.flush_all()
    total_t = time.clock() - total_t

    if output_buffer.num_write_errors > 0 or output_buffer.num_file_errors > 0:
        print >>sys.stderr, '! WARNING: Sweep finished with errors !'
        print >>sys.stderr, '** {writee} reads not written'.format(
            writee=output_buffer.num_write_errors)
        print >>sys.stderr, '** {filee} errors opening files'.format(
            filee=output_buffer.num_file_errors)

    print >>sys.stderr, 'swept {n_reads} for labels...'.format(
        n_reads=n_labeled + n_mlabeled + n_orphaned)
    print >>sys.stderr, '...with {nc} labeled and {no} orphaned'.format(
        nc=n_labeled, no=n_orphaned)
    print >>sys.stderr, '...and {nmc} multilabeled'.format(nmc=n_mlabeled)

    print >>sys.stderr, '** outputting label number distribution...'
    with open('label_dist.txt', 'wb') as outfp:
        for nc in label_number_dist:
            outfp.write('{nc}\n'.format(nc=nc))

if __name__ == '__main__':
    main()
