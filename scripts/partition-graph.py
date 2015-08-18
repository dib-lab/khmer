#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name, missing-docstring
"""
Partition a graph.

% python scripts/partition-graph.py <base>

This will output many <base>.subset.N.pmap files.

Use '-h' for parameter help.
"""
from __future__ import print_function

import threading
import gc
import os.path
import argparse
import textwrap
import sys

import khmer
from khmer.khmer_args import (add_threading_args, info, sanitize_help,
                              ComboFormatter)
from khmer.kfile import check_input_files

# stdlib queue module was renamed on Python 3
try:
    import queue
except ImportError:
    import Queue as queue

DEFAULT_SUBSET_SIZE = int(1e5)
DEFAULT_N_THREADS = 4


def worker(queue, basename, stop_big_traversals):
    while True:
        try:
            (nodegraph, index, start, stop) = queue.get(False)
        except queue.Empty:
            print('exiting', file=sys.stderr)
            return

        outfile = basename + '.subset.%d.pmap' % (index,)
        if os.path.exists(outfile):
            print('SKIPPING', outfile, ' -- already exists', file=sys.stderr)
            continue

        print('starting:', basename, index, file=sys.stderr)

        # pay attention to stoptags when partitioning; take command line
        # direction on whether or not to exhaustively traverse.
        subset = nodegraph.do_subset_partition(start, stop, True,
                                               stop_big_traversals)

        print('saving:', basename, index, file=sys.stderr)
        nodegraph.save_subset_partitionmap(subset, outfile)
        del subset
        gc.collect()


def get_parser():
    epilog = """\
    The resulting partition maps are saved as ``${basename}.subset.#.pmap``
    files.
    """
    parser = argparse.ArgumentParser(
        description="Partition a sequence graph based upon waypoint "
        "connectivity", epilog=textwrap.dedent(epilog), 
        formatter_class=ComboFormatter)

    parser.add_argument('basename', help="basename of the input k-mer"
                        "nodegraph  + tagset files")
    parser.add_argument('--stoptags', '-S', metavar='filename', default='',
                        help="Use stoptags in this file during partitioning")
    parser.add_argument('--subset-size', '-s', default=DEFAULT_SUBSET_SIZE,
                        type=float, help='Set subset size (usually 1e5-1e6 is '
                        'good)')
    parser.add_argument('--no-big-traverse', action='store_true',
                        default=False, help='Truncate graph joins at big '
                        'traversals')
    parser.add_argument('--version', action='version', version='%(prog)s ' +
                        khmer.__version__)
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Overwrite output file if it exists')
    add_threading_args(parser)
    return parser


def main():
    info('partition-graph.py', ['graph'])
    args = sanitize_help(get_parser()).parse_args()
    basename = args.basename

    filenames = [basename, basename + '.tagset']
    for _ in filenames:
        check_input_files(_, args.force)

    print('--', file=sys.stderr)
    print('SUBSET SIZE', args.subset_size, file=sys.stderr)
    print('N THREADS', args.threads, file=sys.stderr)
    if args.stoptags:
        print('stoptag file:', args.stoptags, file=sys.stderr)
    print('--', file=sys.stderr)

    print('loading nodegraph %s' % basename, file=sys.stderr)
    nodegraph = khmer.load_nodegraph(basename)
    nodegraph.load_tagset(basename + '.tagset')

    # do we want to load stop tags, and do they exist?
    if args.stoptags:
        print('loading stoptags from', args.stoptags, file=sys.stderr)
        nodegraph.load_stop_tags(args.stoptags)

    # do we want to exhaustively traverse the graph?
    stop_big_traversals = args.no_big_traverse
    if stop_big_traversals:
        print('** This script brakes for lumps:',
              ' stop_big_traversals is true.', file=sys.stderr)
    else:
        print('** Traverse all the things:',
              ' stop_big_traversals is false.', file=sys.stderr)

    #
    # now, partition!
    #

    # divide the tags up into subsets
    divvy = nodegraph.divide_tags_into_subsets(int(args.subset_size))
    n_subsets = len(divvy)
    divvy.append(0)

    # build a queue of tasks:
    worker_q = queue.Queue()

    # break up the subsets into a list of worker tasks
    for _ in range(0, n_subsets):
        start = divvy[_]
        end = divvy[_ + 1]
        worker_q.put((nodegraph, _, start, end))

    print('enqueued %d subset tasks' % n_subsets, file=sys.stderr)
    open('%s.info' % basename, 'w').write('%d subsets total\n' % (n_subsets))

    n_threads = args.threads
    if n_subsets < n_threads:
        n_threads = n_subsets

    # start threads!
    print('starting %d threads' % n_threads, file=sys.stderr)
    print('---', file=sys.stderr)

    threads = []
    for _ in range(n_threads):
        cur_thrd = threading.Thread(target=worker, args=(worker_q, basename,
                                                         stop_big_traversals))
        threads.append(cur_thrd)
        cur_thrd.start()

    print('done starting threads', file=sys.stderr)

    # wait for threads
    for _ in threads:
        _.join()

    print('---', file=sys.stderr)
    print('done making subsets! see %s.subset.*.pmap' %
          (basename,), file=sys.stderr)

if __name__ == '__main__':
    main()
