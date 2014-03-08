#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import khmer
import sys
import screed
import threading
import Queue

K = 14
HASHTABLE_SIZE = int(1e9)
THRESHOLD = 100

GROUPSIZE = 100
WORKER_THREADS = 4


def process(inq, outq, ht):
    global worker_count

    while not done or not inq.empty():
        try:
            recordlist = inq.get(True, 1)
        except Queue.Empty:
            continue

        x = []
        for record in recordlist:
            kmer = record['sequence'][:K]
            size = ht.calc_connected_graph_size(kmer, THRESHOLD)
            if size >= THRESHOLD:
                x.append(record)

        outq.put(x)

    worker_count -= 1


def write(outq, outfp):
    global worker_count
    while worker_count > 0 or not outq.empty():
        try:
            recordlist = outq.get(True, 1)
        except Queue.Empty:
            continue

        for record in recordlist:
            outfp.write('>%s\n%s\n' % (record['name'], record['sequence']))


def main():
    global done, worker_count
    done = False
    worker_count = 0

    infile = sys.argv[1]
    outfile = infile + '.graphsize2'

    print 'creating ht'
    ht = khmer.new_hashbits(K, HASHTABLE_SIZE, 1)
    print 'eating fa', infile
    total_reads, n_consumed = ht.consume_fasta(infile)
    outfp = open(outfile, 'w')

    inqueue = Queue.Queue(50)
    outqueue = Queue.Queue(50)

    # worker and writer threads
    for i in range(WORKER_THREADS):
        t = threading.Thread(target=process, args=(inqueue, outqueue, ht))
        worker_count += 1
        t.start()

    threading.Thread(target=write, args=(outqueue, outfp)).start()

    # main thread
    x = []
    i = 0
    for n, record in enumerate(screed.fasta.fasta_iter(open(infile))):
        if n % 10000 == 0:
            print '...', n

        x.append(record)
        i += 1

        if i > GROUPSIZE:
            inqueue.put(x)
            x = []
            i = 0
    inqueue.put(x)

    done = True

if __name__ == '__main__':
    main()
