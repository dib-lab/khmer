#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import khmer, sys
import threading
import Queue
import gc
import os.path
import pickle
import math
import random
import subprocess

import glob
import os

K=32

def calc_ht_size(m, k):
   return int(m / k)

def calc_m(n, p):
   n = float(n)
   p = float(p)

   return int(0 - (n*math.log(p))/(math.log(2)**2))

def opt_ht(m, n):
   m = float(m)
   n = float(n)

   k = (m / n) * math.log(2)

   return int(max(1, round(k)))

p = float(sys.argv[2])
n_kmers = 9969000
m = calc_m(n_kmers, p)
N_HT = opt_ht(m, n_kmers)
HASHTABLE_SIZE=calc_ht_size(m, N_HT)

SUBSET_SIZE = int(1e4)
N_THREADS = 1

ht = None

###

save_ht = False
load_ht = False
save_merged_pmap = True
remove_orig_pmap = False

stop_after_n_subsets = None     # only do this many subsets (None == do all)

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
    
    n = 5

    basename = os.path.basename(filename)

    fd = open("log.txt", "w")
    
    primes = []

    below = khmer.get_n_primes_near_x(N_HT * n, HASHTABLE_SIZE)
    above = khmer.get_n_primes_above_x(N_HT * n, HASHTABLE_SIZE)

    primes = below + above
    random.shuffle(primes)

    for run in range(n):

        print primes[run*N_HT:run*N_HT+N_HT]

        ht = khmer._new_hashbits(K, primes[run*N_HT:run*N_HT+N_HT])
        #ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)
    
        # populate the hash table and tag set
        if not load_ht:
            ht.consume_fasta_and_tag(filename)

            # save to a file (optional)
            if save_ht:
                ht.save(basename + '.ht')
                ht.save_tagset(basename + '.tagset')
            
            # calculate the hashtable occupancy
            print '---'
            print 'hashtable occupancy:', ht.n_occupied() / float(HASHTABLE_SIZE)
            print '---'
        else:
            ht.load(basename + '.ht')
            ht.load_tagset(basename + '.tagset')

        # did we just want to load the ht/tagset?
        if stop_after_n_subsets == 0:
            sys.exit(0)

        #stop_tags = pickle.load(open(sys.argv[2]))

        #for stop_tag in stop_tags:
        #    ht.add_stop_tag(stop_tag)

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
            end = divvy[i+1]
            worker_q.put((ht, i, start, end))

        open('%s.info' % basename, 'w').write('%d subsets total\n' % (n_subsets))

        threads = []
        for th in range(N_THREADS):
            t = threading.Thread(target=worker, args=(worker_q, basename))
            threads.append(t)
            t.start()

        # wait for threads
        for t in threads:
            t.join()

        ###

        del ht
        gc.collect()

        # create a new, empty ht object for merging; K matters, but not
        # hashtable size.
        ht = khmer.new_hashbits(K, 1, 1)

        # load & merge all pmap files
        for i in range(0, n_subsets):
            pmap_file = basename + '.subset.%d.pmap' % (i,)
            ht.merge_subset_from_disk(pmap_file)

        # save merged partitionmap
        if save_merged_pmap:
            ht.save_partitionmap(basename + '.pmap.merged')

        if remove_orig_pmap:
            for i in range(0, n_subsets):
                pmap_file = basename + '.subset.%d.pmap' % (i,)
                os.unlink(pmap_file)

        # output partitions!
        n_partitions = ht.output_partitions(filename, basename + '.part')
        (n_partitions, n_singletons) = ht.count_partitions()
        print n_partitions

        fd.write(str(n_partitions) + "\n")
        #print os.listdir(os.getcwd())

        for file in glob.glob(os.getcwd() + "/*pmap*"):
            os.remove(file)
        for file in glob.glob(os.getcwd() + "/*.info"):
            os.remove(file)
        for file in glob.glob(os.getcwd() + "/*.part"):
            os.remove(file)

    fd.close()

if __name__ == '__main__':
    main(sys.argv[1])
