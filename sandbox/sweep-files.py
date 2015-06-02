#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE. Contact: ctb@msu.edu
#
# pylint: disable=invalid-name,missing-docstring,no-member

"""
Find all reads connected to the given contigs on a per-partition basis.

% sweep-files.py -r <range> --db <fasta/q files> \
--query <fasta/q files separate>
"""
from __future__ import print_function

EPILOG = """
Output will be a collection of fasta/q files, each corresponding to a database
file. Each of these is the subset of sequences in the query files which are 
connected to the sequences in the given database file in the de Bruijn graph.
The --range flag sets the breadth of traversal when searching for matches
with the database sequences; for general use, the default (maximum range)
is recommended.

By default, the script will use no output prefix and put outputs in the current
directory. This behavior can be changed with the --prefix and --outdir flags.

The script also uses a queue to improve IO performance with many files. By
default it buffers 1000 reads at a time; this can be changed with 
--max_queue_size.
"""

import screed
import sys
from collections import defaultdict, deque
import os
import time
import khmer
from khmer.khmer_args import (build_hashbits_args, report_on_config, info)

DEFAULT_OUT_PREF = 'reads'
DEFAULT_RANGE = -1

MIN_HSIZE = 4e7
MIN_KSIZE = 21


def get_parser():
    parser = build_hashbits_args('Takes a partitioned reference file \
                                  and a list of reads, and sorts reads \
                                  by which partition they connect to')
    parser.epilog = EPILOG
    parser.add_argument(
        '-r', '--traversal_range', type=int, dest='traversal_range',
        default=DEFAULT_RANGE, help='depth of breadth-first search to perform\
                                    from each read')
    parser.add_argument('--max_queue_size', type=int, default=1000)
    parser.add_argument('--prefix', dest='output_prefix',
                        default=DEFAULT_OUT_PREF,
                        help='Prefix for sorted read files')
    parser.add_argument('--outdir', dest='outdir', default='',
                        help='output directory; default is location of \
                              fastp file')
    parser.add_argument('--query', dest='query', nargs='+',
                        help='Reads to be swept and sorted')
    parser.add_argument('--db', dest='db', nargs='+',
                        help='Database reads for sweep', required=True)

    return parser


def output_single(r):
    if hasattr(r, 'quality'):
        return "@%s\n%s\n+\n%s\n" % (r.name, r.sequence, r.quality)
    else:
        return ">%s\n%s\n" % (r.name, r.sequence)


'''
Simple deque subclass for flushing to a file when a maximum size
is exceeded.
'''
class IODeque(deque):

    def __init__(self, limit, outfp):
        deque.__init__(self)
        self.outfp = outfp
        self.limit = limit

    def append(self, x):
        deque.append(self, x)
        if len(self) >= self.limit:
            self.clear()

    def clear(self):
        while(len(self)):
            self.outfp.write(output_single(self.popleft()))
        deque.clear(self)


def main():
    #info('sweep-files.py', ['sweep'])
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

    outputs = {}

    # Consume the database files and assign each a unique label in the
    # de Bruin graph; open a file and output queue for each file as well.
    ht = khmer.LabelHash(K, HT_SIZE, N_HT)
    try:
        print('consuming and labeling input sequences...', file=sys.stderr)

        for i, dbfile in enumerate(args.db):

            name = args.output_prefix + os.path.basename(dbfile)
            outfp = open(os.path.join(args.outdir, name) + '.sweep', 'wb')
            outq = IODeque(args.max_queue_size, outfp)
            outputs[i] = outq

            for n, record in enumerate(screed.open(dbfile)):
                if n % 50000 == 0:
                    print('...consumed {n} sequences...'.format(n=n), file=sys.stderr)
                ht.consume_sequence_and_tag_with_labels(record.sequence, i)


    except IOError as e:
        print('!! ERROR: !!', e, file=sys.stderr)
        print('...error setting up outputs. exiting...', file=sys.stderr)

    print('done consuming input sequence. \
                        added {t} tags and {l} labels...' \
                        .format(t=ht.n_tags(), l=ht.n_labels()), file=sys.stderr)

    n_orphaned = 0
    n_labeled = 0
    n_mlabeled = 0

    # Iterate through all the reads and check for the labels with which they
    # intersect. Queue to the corresponding label when found.
    for read_file in args.query:
        print('** sweeping {read_file} for labels...'.format(
            read_file=read_file), file=sys.stderr)
        try:
            read_fp = screed.open(read_file)
        except IOError as error:
            print('!! ERROR: !!', error, file=sys.stderr)
            print('*** Could not open {fn}, skipping...'.format(
                fn=read_file), file=sys.stderr)
        else:
            for n, record in enumerate(read_fp):
                if n % 50000 == 0 and n > 0:
                    print('\tswept {n} reads [{nc} labeled, {no} orphaned]' \
                                        .format(n=n, nc=n_labeled,
                                                no=n_orphaned), file=sys.stderr)
                seq = record.sequence
                try:
                    labels = ht.sweep_label_neighborhood(seq, traversal_range)
                except ValueError as e:
                    # sweep_label_neighborhood throws a ValueError when
                    # len(seq) < K. just catch it and move on.
                    pass
                else:
                    if labels:
                        n_labeled += 1
                        if len(labels) > 1:
                            n_mlabeled += 1
                        for label in labels:
                            outputs[label].append(record)
                    else:
                        n_orphaned += 1

            print('** End of file {fn}...'.format(fn=read_file), file=sys.stderr)
            read_fp.close()

    # gotta output anything left in the buffers at the end!
    print('** End of run...', file=sys.stderr)
    for q in list(outputs.values()):
        q.clear()

    print('swept {n_reads}...'.format(
        n_reads=n_labeled + n_orphaned), file=sys.stderr)
    print('...with {nc} labeled and {no} orphaned'.format(
        nc=n_labeled, no=n_orphaned), file=sys.stderr)
    print('...and {nmc} multilabeled'.format(nmc=n_mlabeled), file=sys.stderr)

if __name__ == '__main__':
    main()
