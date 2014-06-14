#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
# pylint: disable=invalid-name,missing-docstring,no-member

"""
Find all reads connected to the given contigs on a per-partition basis.

% sweep-reads-buffered.py -r <range> <contigs fastp> \
<reads1> <reads2> ... <readsN>
"""

EPILOG = """
Output will be a collection of files corresponding to the partitions;
each partition gets a file (prefixed with the output prefix option),
which means this could output many tens or hundreds of thousands of files.
Users should plan accordingly.

This script is very lenient on IO errors, due to the large number of file
operations needed. Thus, errors opening a file for buffer flush or writing
a read to a file will not crash the program; instead, if there were errors,
the user will be warned at the end of execution. Errors with opening read files
are also handled -- we move on to the next read file if there is an error
opening.
"""

import screed
import sys
from collections import defaultdict
import os
import time
import khmer
from khmer.khmer_args import (build_hashbits_args, report_on_config, info)
from khmer.file import (check_file_status, check_valid_file_exists,
                        check_space)

DEFAULT_NUM_BUFFERS = 50000
DEFAULT_MAX_READS = 1000000
DEFAULT_BUFFER_SIZE = 10
DEFAULT_OUT_PREF = 'reads'
DEFAULT_RANGE = -1

MIN_HSIZE = 4e7
MIN_KSIZE = 21


def fmt_fasta(name, seq, labels=[]):
    return '>{name}\t{labels}\n{seq}\n'.format(
        name=name, labels='\t'.join([str(l) for l in labels]), seq=seq)


class ReadBuffer(object):

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


class ReadBufferManager(object):

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
        fn = '{prefix}_{buffer_id}.fa'.format(prefix=self.output_pref,
                                              buffer_id=buf_id)
        fpath = os.path.join(self.outdir, fn)
        buf = self.buffers[buf_id]
        try:
            outfp = open(fpath, 'a')
        except IOError as _:
            print >>sys.stderr, '!! ERROR: {_} !!'.format(_=_)
            print >>sys.stderr, '*** Failed to open {fn} for \
                                buffer flush'.format(fn=fpath)
            self.num_file_errors += 1
        else:
            outfp.write(buf.flush())
            outfp.close()
        finally:
            self.cur_reads -= len(buf)
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


def get_parser():
    parser = build_hashbits_args('Takes a partitioned reference file \
                                  and a list of reads, and sorts reads \
                                  by which partition they connect to')
    parser.epilog = EPILOG
    parser.add_argument(
        '-r', '--traversal_range', type=int, dest='traversal_range',
        default=DEFAULT_RANGE, help='depth of breadth-first search to perform\
                                    from each read')
    parser.add_argument('-b', '--buffer_size', dest='max_reads', type=int,
                        default=DEFAULT_MAX_READS,
                        help='Max total reads to buffer before flushing')
    parser.add_argument('-l', '--buffer_length', dest='buffer_size', type=int,
                        default=DEFAULT_BUFFER_SIZE,
                        help='Max length of an individual label buffer \
                              before flushing')
    parser.add_argument('--prefix', dest='output_prefix',
                        default=DEFAULT_OUT_PREF,
                        help='Prefix for sorted read files')
    parser.add_argument('--outdir', dest='outdir',
                        help='output directory; default is location of \
                              fastp file')
    parser.add_argument('-m', '--max_buffers', dest='max_buffers', type=int,
                        default=DEFAULT_NUM_BUFFERS,
                        help='Max individual label buffers before flushing')
    labeling = parser.add_mutually_exclusive_group(required=True)
    labeling.add_argument('--label-by-pid', dest='label_by_pid',
                          action='store_true', help='separate reads by\
                        referece partition id')
    labeling.add_argument('--label-by-seq', dest='label_by_seq',
                          action='store_true', help='separate reads by\
                        reference sequence')
    labeling.add_argument('--label-by-group', dest='group_size', type=int,
                          help='separate reads by arbitrary sized groups\
                        of reference sequences')
    parser.add_argument(dest='input_fastp', help='Reference fasta or fastp')
    parser.add_argument('input_files', nargs='+',
                        help='Reads to be swept and sorted')

    return parser


def main():
    info('sweep-reads-buffered.py', ['sweep'])
    parser = get_parser()
    args = parser.parse_args()

    if args.min_tablesize < MIN_HSIZE:
        args.min_tablesize = MIN_HSIZE
    if args.ksize < MIN_KSIZE:
        args.ksize = MIN_KSIZE

    report_on_config(args, hashtype='hashbits')

    K = args.ksize
    HT_SIZE = args.min_tablesize
    N_HT = args.n_tables

    traversal_range = args.traversal_range
    input_fastp = args.input_fastp

    if not args.outdir:
        outdir = os.path.dirname(input_fastp)
    else:
        outdir = args.outdir

    max_buffers = args.max_buffers
    output_pref = args.output_prefix
    buf_size = args.buffer_size
    max_reads = args.max_reads

    check_file_status(args.input_fastp)
    check_valid_file_exists(args.input_files)
    all_input_files = [input_fastp]
    all_input_files.extend(args.input_files)
    # Check disk space availability
    check_space(all_input_files)

    output_buffer = ReadBufferManager(
        max_buffers, max_reads, buf_size, output_pref, outdir)

    # consume the partitioned fasta with which to label the graph
    ht = khmer.LabelHash(K, HT_SIZE, N_HT)
    try:
        print >>sys.stderr, 'consuming input sequences...'
        if args.label_by_pid:
            print >>sys.stderr, '...labeling by partition id (pid)'
            ht.consume_partitioned_fasta_and_tag_with_labels(input_fastp)
        elif args.label_by_seq:
            print >>sys.stderr, '...labeling by sequence'
            for n, record in enumerate(screed.open(input_fastp)):
                if n % 50000 == 0:
                    print >>sys.stderr, \
                        '...consumed {n} sequences...'.format(n=n)
                ht.consume_sequence_and_tag_with_labels(record.sequence, n)
        else:
            print >>sys.stderr, \
                '...labeling to create groups of size {s}'.format(
                    s=args.group_size)
            label = -1
            g = 0
            try:
                outfn = os.path.join(args.outdir, 
                    '{pref}_base_{g}.fa'.format(pref=output_pref, g=g) 
                outfp = open(outfn, 'wb')
                for n, record in enumerate(screed.open(input_fastp)):
                    if n % args.group_size == 0:
                        label += 1
                        if label > g:
                            g = label
                            outfp = open('{pref}_base_{g}.fa'.format(
                                pref=output_pref, g=g), 'wb')
                    if n % 50000 == 0:
                        print >>sys.stderr, \
                            '...consumed {n} sequences...'.format(n=n)
                    ht.consume_sequence_and_tag_with_labels(record.sequence,
                                                            label)
                    outfp.write('>{name}\n{seq}\n'.format(name=record.name,
                                                          seq=record.sequence))

            except IOError as e:
                print >>sys.stderr, '!! ERROR !!', e
                print >>sys.stderr, '...error splitting input. exiting...'

    except IOError as e:
        print >>sys.stderr, '!! ERROR: !!', e
        print >>sys.stderr, '...error consuming \
                            {i}. exiting...'.format(i=input_fastp)

    print >>sys.stderr, 'done consuming input sequence. \
                        added {t} tags and {l} \
                        labels...'.format(t=ht.n_tags(), l=ht.n_labels())

    label_dict = defaultdict(int)
    label_number_dist = []

    n_orphaned = 0
    n_labeled = 0
    n_mlabeled = 0

    total_t = time.clock()
    start_t = time.clock()
    for read_file in args.input_files:
        print >>sys.stderr, '** sweeping {read_file} for labels...'.format(
            read_file=read_file)
        file_t = 0.0
        try:
            read_fp = screed.open(read_file)
        except IOError as error:
            print >>sys.stderr, '!! ERROR: !!', error
            print >>sys.stderr, '*** Could not open {fn}, skipping...'.format(
                fn=read_file)
        else:
            for _, record in enumerate(read_fp):
                if _ % 50000 == 0:
                    end_t = time.clock()
                    batch_t = end_t - start_t
                    file_t += batch_t
                    print >>sys.stderr, '\tswept {n} reads [{nc} labeled, \
                                         {no} orphaned] \
                                        ** {sec}s ({sect}s total)' \
                                        .format(n=_, nc=n_labeled,
                                                no=n_orphaned,
                                                sec=batch_t, sect=file_t)
                    start_t = time.clock()
                seq = record.sequence
                name = record.name
                try:
                    labels = ht.sweep_label_neighborhood(seq, traversal_range)
                except ValueError as e:
                    pass
                else:
                    seq_str = fmt_fasta(name, seq, labels)
                    label_number_dist.append(len(labels))
                    if labels:
                        n_labeled += 1
                        if len(labels) > 1:
                            output_buffer.queue(seq_str, 'multi')
                            n_mlabeled += 1
                            label_dict['multi'] += 1
                        else:
                            output_buffer.queue(seq_str, labels[0])
                            label_dict[labels[0]] += 1
                    else:
                        n_orphaned += 1
                        output_buffer.queue(seq_str, 'orphaned')
                        label_dict['orphaned'] += 1
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
        n_reads=n_labeled + n_orphaned)
    print >>sys.stderr, '...with {nc} labeled and {no} orphaned'.format(
        nc=n_labeled, no=n_orphaned)
    print >>sys.stderr, '...and {nmc} multilabeled'.format(nmc=n_mlabeled)

    print >>sys.stderr, '** outputting label number distribution...'
    fn = os.path.join(outdir, '{pref}.dist.txt'.format(pref=output_pref))
    with open(fn, 'wb') as outfp:
        for nc in label_number_dist:
            outfp.write('{nc}\n'.format(nc=nc))

    fn = os.path.join(outdir, '{pref}.counts.csv'.format(pref=output_pref))
    print >>sys.stderr, '** outputting label read counts...'
    with open(fn, 'wb') as outfp:
        for k in label_dict:
            outfp.write('{l},{c}\n'.format(l=k, c=label_dict[k]))

if __name__ == '__main__':
    main()
