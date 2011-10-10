# /usr/bin/env python
"""
Partition a graph.

% python scripts/partition-graph.py <base>

Use '-h' for parameter help.
"""

import khmer, sys
import threading
import Queue
import gc
import os.path
import argparse

DEFAULT_SUBSET_SIZE = int(1e5)
DEFAULT_N_THREADS = 4

ht = None

###

save_ht = True
load_ht = False
save_merged_pmap = True
remove_orig_pmap = False

assert not (save_ht and load_ht)         # incompatible

if not save_merged_pmap and remove_orig_pmap:
    print '** warning, all the pmap files are going away! no permanent record!'
    print ''

###

def worker(q, basename):
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

        # pay attention to stoptags when partitioning, note
        subset = ht.do_subset_partition(start, stop, True)

        print 'saving:', basename, n
        ht.save_subset_partitionmap(subset, outfile)
        del subset
        gc.collect()

def main():
    parser = argparse.ArgumentParser(description="Partition a graph.")

    parser.add_argument('basename')
    parser.add_argument('--stoptags', '-S', dest='stoptags', default='',
                        help="Use stoptags in this file during partitioning")
    parser.add_argument('--no-merge', '-M', dest='merge_subsets', default=True,
                        action='store_false', help='Do not merge subsets.')
    parser.add_argument('--keep-subsets', dest='remove_subsets',
                        default=True, action='store_false',
                        help='Keep individual subsets')
    parser.add_argument('--subset-size', '-s', default=DEFAULT_SUBSET_SIZE,
                        dest='subset_size', type=float,
                        help='Set subset size (usually 1e5-1e6 is good)')

    parser.add_argument('--threads', '-T', dest='n_threads',
                        default=DEFAULT_N_THREADS,
                        help='Number of simultaneous threads to execute')

    args = parser.parse_args()
    basename = args.basename

    print '--'
    print 'SUBSET SIZE', args.subset_size
    print 'N THREADS', args.n_threads
    print 'merge subsets?', args.merge_subsets
    print 'remove subsets after merging?', args.remove_subsets
    if args.stoptags:
        print 'stoptag file:', args.stoptags
    print '--'

    print 'loading ht %s.ht' % basename
    ht = khmer.load_hashbits(basename + '.ht')
    ht.load_tagset(basename + '.tagset')

    # retrieve K
    K = ht.ksize()

    # do we want to load stop tags, and do they exist? @CTB
    stoptags_file = basename + '.stoptags'
    if args.stoptags:
        print 'loading stoptags from', args.stoptags
        ht.load_stop_tags(stoptags_file)

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

    n_threads = args.n_threads
    if n_subsets < n_threads:
        n_threads = n_subsets

    # start threads!
    print 'starting %d threads' % n_threads
    print '---'

    threads = []
    for n in range(n_threads):
        t = threading.Thread(target=worker, args=(worker_q, basename))
        threads.append(t)
        t.start()

    print 'done starting threads'

    # wait for threads
    for t in threads:
        t.join()

    print '---'
    print 'done making subsets! see %s.subset.*.pmap' % (basename,)

    ###

    # load & merge all pmap files
    if args.merge_subsets:
        print 'erasing old ht, creating new for merge'
        del ht
        gc.collect()

        # create a new, empty ht object for merging; K matters, but not
        # hashtable size.
        ht = khmer.new_hashbits(K, 1, 1)

        for i in range(0, n_subsets):
            pmap_file = basename + '.subset.%d.pmap' % (i,)
            print 'merging', pmap_file
            ht.merge_subset_from_disk(pmap_file)

        # save merged partitionmap
        print 'saving merged pmap to %s.pmap.merged' % basename 
        ht.save_partitionmap(basename + '.pmap.merged')

        if args.remove_subsets:
            print 'removing subset pmap files'
            for i in range(0, n_subsets):
                pmap_file = basename + '.subset.%d.pmap' % (i,)
                os.unlink(pmap_file)

        print '---'
        print 'see final partitioning in', basename + '.pmap.merged'

if __name__ == '__main__':
    main()
