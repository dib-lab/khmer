#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import khmer
import sys
import threading
import Queue
import gc
import os.path

K = 32
HASHTABLE_SIZE = int(8e8)
N_HT = 4

SUBSET_SIZE = int(1e4)
N_THREADS = 8

ht = None

###

save_ht = True
load_ht = False
save_merged_pmap = True
remove_orig_pmap = False

stop_after_n_subsets = None     # only do this many subsets (None == do all)
load_stoptags_if_exist = True   # load stoptags, if a .stoptags file exists.

assert not (save_ht and load_ht)         # incompatible
if stop_after_n_subsets == 0:
    assert save_ht

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
            print 'SKIPPING', basename, ' -- already exists'
            continue

        print 'starting:', basename, n

        # pay attention to stoptags when partitioning, note
        subset = ht.do_subset_partition(start, stop, True)

        print 'saving:', basename, n
        ht.save_subset_partitionmap(subset, outfile)
        del subset
        gc.collect()


def main(filename):
    global ht

    basename = os.path.basename(filename)

    print 'input file to partition: %s' % filename
    print '-- settings:'
    print 'K', K
    print 'HASHTABLE SIZE %g' % HASHTABLE_SIZE
    print 'N HASHTABLES %d' % N_HT
    print 'SUBSET SIZE', SUBSET_SIZE
    print 'N THREADS', N_THREADS
    print '--'

    ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)

    # populate the hash table and tag set
    if not load_ht:
        print 'reading sequences and loading tagset from %s...' % (filename,)
        ht.consume_fasta_and_tag(filename)

        # save to a file (optional)
        if save_ht:
            print 'saving...'
            ht.save(basename + '.ht')
            print 'saving tagset...'
            ht.save_tagset(basename + '.tagset')

        # calculate the hashtable occupancy
        print '---'
        print 'hashtable occupancy:', ht.n_occupied() / float(HASHTABLE_SIZE)
        print '---'
    else:
        print 'loading ht %s.ht' % basename
        ht.load(basename + '.ht')
        print 'loading tagset %s.tagset...' % basename
        ht.load_tagset(basename + '.tagset')

    # did we just want to load the ht/tagset?
    if stop_after_n_subsets == 0:
        sys.exit(0)

    # do we want to load stop tags, and do they exist?
    stoptags_file = basename + '.stoptags'
    if load_stoptags_if_exist and os.path.exists(stoptags_file):
        print 'loading stoptags from', stoptags_file
        ht.load_stop_tags(stoptags_file)

    #
    # now, partition!
    #

    # divide the tags up into subsets
    divvy = ht.divide_tags_into_subsets(SUBSET_SIZE)
    n_subsets = len(divvy)
    divvy.append(0)

    # build a queue of tasks:
    worker_q = Queue.Queue()

    for i in range(0, n_subsets):
        if stop_after_n_subsets is not None and i >= stop_after_n_subsets:
            break

        start = divvy[i]
        end = divvy[i + 1]
        worker_q.put((ht, i, start, end))

    print 'enqueued %d subset tasks' % n_subsets
    open('%s.info' % basename, 'w').write('%d subsets total\n' % (n_subsets))

    threads = []
    for n in range(N_THREADS):
        t = threading.Thread(target=worker, args=(worker_q, basename))
        threads.append(t)
        t.start()

    print 'started threads'

    # wait for threads
    for t in threads:
        t.join()

    print 'done making subsets! see %s.subset.*.pmap' % (basename,)

    ###

    print 'erasing old ht, creating new'
    del ht
    gc.collect()

    # create a new, empty ht object for merging; K matters, but not
    # hashtable size.
    ht = khmer.new_hashbits(K, 1, 1)

    # load & merge all pmap files
    for i in range(0, n_subsets):
        pmap_file = basename + '.subset.%d.pmap' % (i,)
        print 'loading', pmap_file
        ht.merge_subset_from_disk(pmap_file)

    # save merged partitionmap
    if save_merged_pmap:
        print 'saving merged pmap to %s.pmap.merged' % basename
        ht.save_partitionmap(basename + '.pmap.merged')

    if remove_orig_pmap:
        print 'removing subset pmap files'
        for i in range(0, n_subsets):
            pmap_file = basename + '.subset.%d.pmap' % (i,)
            os.unlink(pmap_file)

    # output partitions!
    n_partitions = ht.output_partitions(filename, basename + '.part')
    (n_partitions, n_singletons) = ht.count_partitions()
    print 'output partitions:', n_partitions
    print 'pmap partitions:', n_partitions
    print 'singletons:', n_singletons


if __name__ == '__main__':
    main(sys.argv[1])
