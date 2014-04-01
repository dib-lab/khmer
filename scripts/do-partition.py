#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=missing-docstring,invalid-name
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
import textwrap
from khmer.khmer_args import build_hashbits_args
from khmer.khmer_args import report_on_config
import glob
from khmer.file import check_file_status, check_space

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
    def __debug_vm_usage(msg):  # pylint: disable=unused-argument
        pass


def worker(queue, basename, stop_big_traversals):
    while 1:
        try:
            (htable, index, start, stop) = queue.get(False)
        except Queue.Empty:
            print 'exiting'
            return

        outfile = basename + '.subset.%d.pmap' % (index,)
        if os.path.exists(outfile):
            print 'SKIPPING', outfile, ' -- already exists'
            continue

        print 'starting:', basename, index

        # pay attention to stoptags when partitioning; take command line
        # direction on whether or not to exhaustively traverse.
        subset = htable.do_subset_partition(start, stop, True,
                                            stop_big_traversals)

        print 'saving:', basename, index
        htable.save_subset_partitionmap(subset, outfile)
        del subset
        gc.collect()


def get_parser():
    epilog = """
    Load in a set of sequences, partition them, merge the partitions, and
    annotate the original sequences files with the partition information.

    This script combines the functionality of :program:`load-graph.py`,
    :program:`partition-graph.py`, :program:`merge-partitions.py`, and
    :program:`annotate-partitions.py` into one script. This is convenient
    but should probably not be used for large data sets, because
    :program:`do-partition.py` doesn't provide save/resume functionality.
    """
    parser = build_hashbits_args(
        descr='Load, partition, and annotate FAST[AQ] sequences',
        epilog=textwrap.dedent(epilog))
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
    parser.add_argument('graphbase', help="base name for output files")
    parser.add_argument('input_filenames', nargs='+')
    return parser


# pylint: disable=too-many-branches
def main():  # pylint: disable=too-many-locals,too-many-statements
    args = get_parser().parse_args()

    report_on_config(args, hashtype='hashbits')

    for infile in args.input_filenames:
        check_file_status(infile)

    check_space(args.input_filenames)

    print 'Saving k-mer presence table to %s' % args.graphbase
    print 'Loading kmers from sequences in %s' % repr(args.input_filenames)

    print '--'
    print 'SUBSET SIZE', args.subset_size
    print 'N THREADS', args.n_threads
    print '--'

    # load-graph

    print 'making k-mer presence table'
    htable = khmer.new_hashbits(args.ksize, args.min_tablesize, args.n_tables)

    for _, filename in enumerate(args.input_filenames):
        print 'consuming input', filename
        htable.consume_fasta_and_tag(filename)

    fp_rate = khmer.calc_expected_collisions(htable)
    print 'fp rate estimated to be %1.3f' % fp_rate
    if fp_rate > 0.15:          # 0.18 is ACTUAL MAX. Do not change.
        print >> sys.stderr, "**"
        print >> sys.stderr, ("** ERROR: the graph structure is too small for"
                              " this data set.  Increase k-mer presence table "
                              "size/num of tables.")
        print >> sys.stderr, "**"
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
    divvy = htable.divide_tags_into_subsets(int(args.subset_size))
    n_subsets = len(divvy)
    divvy.append(0)

    # build a queue of tasks:
    worker_q = Queue.Queue()

    # break up the subsets into a list of worker tasks
    for _ in range(0, n_subsets):
        start = divvy[_]
        end = divvy[_ + 1]
        worker_q.put((htable, _, start, end))

    print 'enqueued %d subset tasks' % n_subsets
    open('%s.info' % args.graphbase, 'w').write('%d subsets total\n'
                                                % (n_subsets))

    if n_subsets < args.n_threads:
        args.n_threads = n_subsets

    # start threads!
    print 'starting %d threads' % args.n_threads
    print '---'

    threads = []
    for _ in range(args.n_threads):
        cur_thread = threading.Thread(target=worker,
                                      args=(worker_q, args.graphbase,
                                            stop_big_traversals))
        threads.append(cur_thread)
        cur_thread.start()

    print 'done starting threads'

    # wait for threads
    for _ in threads:
        _.join()

    print '---'
    print 'done making subsets! see %s.subset.*.pmap' % (args.graphbase,)

    # merge-partitions

    pmap_files = glob.glob(args.graphbase + '.subset.*.pmap')

    print 'loading %d pmap files (first one: %s)' % (len(pmap_files),
                                                     pmap_files[0])

    htable = khmer.new_hashbits(args.ksize, 1, 1)

    for pmap_file in pmap_files:
        print 'merging', pmap_file
        htable.merge_subset_from_disk(pmap_file)

    if args.remove_subsets:
        print 'removing pmap files'
        for pmap_file in pmap_files:
            os.unlink(pmap_file)

    # annotate-partitions

    for infile in args.input_filenames:
        print 'outputting partitions for', infile
        outfile = os.path.basename(infile) + '.part'
        part_count = htable.output_partitions(infile, outfile)
        print 'output %d partitions for %s' % (part_count, infile)
        print 'partitions are in', outfile

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
