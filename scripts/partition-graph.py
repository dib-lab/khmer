#! /usr/bin/env python
"""
Partition a graph.

% python scripts/partition-graph.py <base>

This will output many <base>.subset.N.pmap files.

Use '-h' for parameter help.
"""

import sys
import threading
import Queue
import gc
import os.path
import argparse

# Debugging Support
import re
import platform
if "Linux" == platform.system( ):
    def __debug_vm_usage( msg ):
	print "===> DEBUG: " + msg
	for vmstat in re.findall( r".*Vm.*", file( "/proc/self/status" ).read( ) ):
	    print vmstat
else:
    def __debug_vm_usage( msg ): pass

import khmer

DEFAULT_SUBSET_SIZE = int(1e5)
DEFAULT_N_THREADS = 4

###

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
    parser = argparse.ArgumentParser(description="Partition a graph.")

    parser.add_argument('basename')
    parser.add_argument('--stoptags', '-S', dest='stoptags', default='',
                        help="Use stoptags in this file during partitioning")
    parser.add_argument('--subset-size', '-s', default=DEFAULT_SUBSET_SIZE,
                        dest='subset_size', type=float,
                        help='Set subset size (usually 1e5-1e6 is good)')

    parser.add_argument('--no-big-traverse', dest='no_big_traverse',
                        action='store_true', default=False,
                        help='Truncate graph joins at big traversals')

    parser.add_argument('--threads', '-T', dest='n_threads',
                        default=DEFAULT_N_THREADS,
                        help='Number of simultaneous threads to execute')

    args = parser.parse_args()
    basename = args.basename

    print '--'
    print 'SUBSET SIZE', args.subset_size
    print 'N THREADS', args.n_threads
    if args.stoptags:
        print 'stoptag file:', args.stoptags
    print '--'

    print 'loading ht %s.ht' % basename
    ht = khmer.load_hashbits(basename + '.ht')
    ht.load_tagset(basename + '.tagset')

    # retrieve K
    K = ht.ksize()

    # do we want to load stop tags, and do they exist?
    if args.stoptags:
        print 'loading stoptags from', args.stoptags
        ht.load_stop_tags(args.stoptags)

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
        end = divvy[i+1]
        worker_q.put((ht, i, start, end))

    print 'enqueued %d subset tasks' % n_subsets
    open('%s.info' % basename, 'w').write('%d subsets total\n' % (n_subsets))

    n_threads = int(args.n_threads)
    if n_subsets < n_threads:
        n_threads = n_subsets

    # start threads!
    print 'starting %d threads' % n_threads
    print '---'

    threads = []
    for n in range(n_threads):
        t = threading.Thread(target=worker, args=(worker_q, basename,
                                                  stop_big_traversals))
        threads.append(t)
        t.start()

    print 'done starting threads'

    # wait for threads
    for t in threads:
        t.join()

    print '---'
    print 'done making subsets! see %s.subset.*.pmap' % (basename,)

if __name__ == '__main__':
    main()
