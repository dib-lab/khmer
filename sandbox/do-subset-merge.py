#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
import sys
import khmer
import gc
import Queue
import threading
import glob
import os.path
import shutil

K = 32

###

lock = threading.Lock()
output_n = 0

###


def pull_pair(q):
    while q.qsize() >= 2:
        try:
            ht, filename1 = q.get(False)
        except Queue.Empty:
            print 'exiting'
            return

        try:
            ht, filename2 = q.get(False)
        except Queue.Empty:
            print 'pushing, then exiting'
            q.put((ht, filename1))
            return

        merge_file = merge(filename1, filename2, ht)
        # q.put((ht, merge_file))   # @CTB


def merge(filename1, filename2, ht):
    global output_n, lock

    lock.acquire()
    next_n = output_n
    output_n += 1
    lock.release()

    merge_filename = os.path.join(dir2, '%s.merge.%d' % (dir2, next_n))
    print 'merge: %s = %s + %s' % (merge_filename, filename1, filename2)
    subset1 = ht.load_subset_partitionmap(filename1)
    ht.merge2_subset_from_disk(subset1, filename2)
    ht.save_subset_partitionmap(subset1, merge_filename + '.pmap')

    return merge_filename + '.pmap'


def main(dir1, dir2, n_threads):
    # detect all of the relevant partitionmap files
    subset_filenames = glob.glob(os.path.join(dir1, '*.pmap'))

    # create empty hashtable structure
    ht = khmer.new_hashbits(K, 1, 1)

    # put jobs on queue
    merge_queue = Queue.Queue()
    for filename in subset_filenames:
        merge_queue.put((ht, filename))

    print 'starting threads'

    threads = []
    for n in range(n_threads):
        t = threading.Thread(target=pull_pair, args=(merge_queue,))
        threads.append(t)
        t.start()

    # wait for threads
    for t in threads:
        t.join()

    # done!

    if merge_queue.qsize() == 1:
        ht, merge_file = merge_queue.get()
        print 'copying', merge_file
        shutil.copy(merge_file, os.path.join(dir2,
                                             os.path.basename(merge_file)))

    assert merge_queue.qsize() == 0

if __name__ == '__main__':
    n_threads = int(sys.argv[1])
    dir1 = sys.argv[2]
    dir2 = sys.argv[3]

    main(dir1, dir2, n_threads)
