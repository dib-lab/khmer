# !/usr/bin/env python2
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
from collections import defaultdict, deque
import os
import time
import khmer
from khmer.khmer_args import (build_hashbits_args, report_on_config, info)
from khmer.file import (check_file_status, check_valid_file_exists,
                        check_space)

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
                        default='',
                        help='Prefix for sorted read files')
    parser.add_argument('--outdir', dest='outdir', default='',
                        help='output directory; default is location of \
                              fastp file')
    parser.add_argument('--query', dest='query', nargs='+',
                        help='Reads to be swept and sorted')
    parser.add_argument('--db', dest='db', nargs='+',
                        help='Database reads for sweep')

    return parser


def output_single(r):
    if hasattr(r, 'accuracy'):
        return "@%s\n%s\n+\n%s\n" % (r.name, r.sequence, r.accuracy)
    else:
        return ">%s\n%s\n" % (r.name, r.sequence)


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
            self.outfp.write(output_single(self.pop()))
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

    # consume the partitioned fasta with which to label the graph
    ht = khmer.LabelHash(K, HT_SIZE, N_HT)
    try:
        print >>sys.stderr, 'consuming and labeling input sequences...'

        for i, dbfile in enumerate(args.db):

            name = args.output_prefix + os.path.basename(dbfile)
            outfp = open(os.path.join(args.outdir, name) + '.sweep', 'wb')
            outq = IODeque(args.max_queue_size, outfp)
            outputs[i] = outq

            for n, record in enumerate(screed.open(dbfile)):
                if n % 50000 == 0:
                    print >>sys.stderr, \
                        '...consumed {n} sequences...'.format(n=n)
                ht.consume_sequence_and_tag_with_labels(record.sequence, i)


    except IOError as e:
        print >>sys.stderr, '!! ERROR: !!', e
        print >>sys.stderr, '...error setting up outputs. exiting...'

    print >>sys.stderr, 'done consuming input sequence. \
                        added {t} tags and {l} labels...' \
                        .format(t=ht.n_tags(), l=ht.n_labels())

    n_orphaned = 0
    n_labeled = 0
    n_mlabeled = 0

    for read_file in args.query:
        print >>sys.stderr, '** sweeping {read_file} for labels...'.format(
            read_file=read_file)
        try:
            read_fp = screed.open(read_file)
        except IOError as error:
            print >>sys.stderr, '!! ERROR: !!', error
            print >>sys.stderr, '*** Could not open {fn}, skipping...'.format(
                fn=read_file)
        else:
            for n, record in enumerate(read_fp):
                if n % 50000 == 0 and n > 0:
                    print >>sys.stderr, \
                        '\tswept {n} reads [{nc} labeled, {no} orphaned]' \
                                        .format(n=n, nc=n_labeled,
                                                no=n_orphaned)
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

            print >>sys.stderr, '** End of file {fn}...'.format(fn=read_file)
            read_fp.close()

    # gotta output anything left in the buffers at the end!
    print >>sys.stderr, '** End of run...'
    for q in outputs.values():
        q.clear()

    print >>sys.stderr, 'swept {n_reads}...'.format(
        n_reads=n_labeled + n_orphaned)
    print >>sys.stderr, '...with {nc} labeled and {no} orphaned'.format(
        nc=n_labeled, no=n_orphaned)
    print >>sys.stderr, '...and {nmc} multilabeled'.format(nmc=n_mlabeled)

if __name__ == '__main__':
    main()
