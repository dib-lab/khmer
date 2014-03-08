#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
"""
Do all the partition steps in one script.

% do-partition.py <graphname> <reads1> [ <reads2> ... ]

Use '-h' for parameter help.
"""

import khmer
import sys
import threading
import Queue
import gc
import os.path
import os
from khmer.khmer_args import build_hashbits_args
from khmer.khmer_args import report_on_config
import glob
from khmer.file_api import check_file_status, check_space

DEFAULT_SUBSET_SIZE = int(1e5)
DEFAULT_N_THREADS = 4
DEFAULT_K = 32

# Debugging Support
import re
import platform
if "Linux" == platform.system():
    def __debug_vm_usage(msg):
        print "===> DEBUG: " + msg
        for vmstat in re.findall(r".*Vm.*", file("/proc/self/status").read()):
            print vmstat
else:
    def __debug_vm_usage(msg):
        pass


def worker(q, basename, stop_big_traversals):
    while 1:
        try:
            (ht, n, start, stop) = q.get(False)
        except Queue.Empty:
            print 'exiting'
            return

        outfile = basename + '.subset.%d.pmap' % (n,)
        if os.path.exists(outfile):
            print 'SKIPPING', outfile, ' -- already exists'
            continue

        print 'starting:', basename, n

        # pay attention to stoptags when partitioning; take command line
        # direction on whether or not to exhaustively traverse.
        subset = ht.do_subset_partition(start, stop, True, stop_big_traversals)

        print 'saving:', basename, n
        ht.save_subset_partitionmap(subset, outfile)
        del subset
        gc.collect()


def main():
    parser = build_hashbits_args()
    parser.add_argument('--subset-size', '-s', default=DEFAULT_SUBSET_SIZE,
                        dest='subset_size', type=float,
                        help='Set subset size (usually 1e5-1e6 is good)')

    parser.add_argument('--no-big-traverse', dest='no_big_traverse',
                        action='store_true', default=False,
                        help='Truncate graph joins at big traversals')

    parser.add_argument('--threads', '-T', dest='n_threads',
                        default=DEFAULT_N_THREADS,
                        help='Number of simultaneous threads to execute')

    parser.add_argument('--keep-subsets', dest='remove_subsets',
                        default=True, action='store_false',
                        help='Keep individual subsets (default: False)')

    parser.add_argument('graphbase')
    parser.add_argument('input_filenames', nargs='+')

    args = parser.parse_args()

    report_on_config(args, hashtype='hashbits')

    K = args.ksize
    HT_SIZE = args.min_hashsize
    N_HT = args.n_hashes

    base = args.graphbase
    filenames = args.input_filenames

    for infile in filenames:
        check_file_status(infile)

    check_space(filenames)

    print 'Saving hashtable to %s' % base
    print 'Loading kmers from sequences in %s' % repr(filenames)

    print '--'
    print 'SUBSET SIZE', args.subset_size
    print 'N THREADS', args.n_threads
    print '--'

    # load-graph

    print 'making hashtable'
    ht = khmer.new_hashbits(K, HT_SIZE, N_HT)

    for n, filename in enumerate(filenames):
        print 'consuming input', filename
        ht.consume_fasta_and_tag(filename)

    fp_rate = khmer.calc_expected_collisions(ht)
    print 'fp rate estimated to be %1.3f' % fp_rate
    if fp_rate > 0.15:          # 0.18 is ACTUAL MAX. Do not change.
        print >>sys.stderr, "**"
        print >>sys.stderr, "** ERROR: the graph structure is too small for"
        print >>sys.stderr, "** this data set.  Increase hashsize/num ht."
        print >>sys.stderr, "**"
        sys.exit(1)

    # partition-graph

    # do we want to exhaustively traverse the graph?
    stop_big_traversals = args.no_big_traverse
    if stop_big_traversals:
        print '** This script brakes for lumps: stop_big_traversals is true.'
    else:
        print '** Traverse all the things: stop_big_traversals is false.'

    #
    # now, partition!
    #

    # divide the tags up into subsets
    divvy = ht.divide_tags_into_subsets(int(args.subset_size))
    n_subsets = len(divvy)
    divvy.append(0)

    # build a queue of tasks:
    worker_q = Queue.Queue()

    # break up the subsets into a list of worker tasks
    for i in range(0, n_subsets):
        start = divvy[i]
        end = divvy[i + 1]
        worker_q.put((ht, i, start, end))

    print 'enqueued %d subset tasks' % n_subsets
    open('%s.info' % base, 'w').write('%d subsets total\n' % (n_subsets))

    n_threads = int(args.n_threads)
    if n_subsets < n_threads:
        n_threads = n_subsets

    # start threads!
    print 'starting %d threads' % n_threads
    print '---'

    threads = []
    for n in range(n_threads):
        t = threading.Thread(target=worker, args=(worker_q, base,
                                                  stop_big_traversals))
        threads.append(t)
        t.start()

    print 'done starting threads'

    # wait for threads
    for t in threads:
        t.join()

    print '---'
    print 'done making subsets! see %s.subset.*.pmap' % (base,)

    # merge-partitions

    pmap_files = glob.glob(args.graphbase + '.subset.*.pmap')

    print 'loading %d pmap files (first one: %s)' % (len(pmap_files),
                                                     pmap_files[0])

    ht = khmer.new_hashbits(K, 1, 1)

    for pmap_file in pmap_files:
        print 'merging', pmap_file
        ht.merge_subset_from_disk(pmap_file)

    if args.remove_subsets:
        print 'removing pmap files'
        for pmap_file in pmap_files:
            os.unlink(pmap_file)

    # annotate-partitions

    for infile in args.input_filenames:
        print 'outputting partitions for', infile
        outfile = os.path.basename(infile) + '.part'
        n = ht.output_partitions(infile, outfile)
        print 'output %d partitions for %s' % (n, infile)
        print 'partitions are in', outfile

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
