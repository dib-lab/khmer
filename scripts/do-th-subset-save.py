import khmer, sys
import threading
import Queue
import gc
import os.path

K=32
HASHTABLE_SIZE=4**15+1

ht = khmer.new_hashtable(K, HASHTABLE_SIZE)

FILENAME=sys.argv[1]
SUBSET_SIZE = 100000
N_THREADS = 4

def worker(q):
    while 1:
        try:
            (ht, n, start, stop) = q.get(False)
        except Queue.Empty:
            print 'exiting'
            return

        outfile = FILENAME + '.subset.%d' % (n,)
        if os.path.exists(outfile + '.pmap'):
            print 'SKIPPING', FILENAME, start, stop
            continue
        
        print 'starting:', FILENAME, n
        subset = ht.do_subset_partition(start, stop)
        print 'saving:', FILENAME, n
        
        outfile = FILENAME + '.subset.%d' % (n,)
        ht.save_subset_partitionmap(subset,
                                    outfile + '.pmap',
                                    outfile + '.surr')
        del subset
        gc.collect()

(total_reads, total_kmers) = ht.consume_fasta_and_tag(FILENAME)
divvy = ht.divide_tags_into_subsets(SUBSET_SIZE)
n_subsets = len(divvy)
divvy.append(0)

print '---'
print 'hashtable occupancy:', ht.n_occupied() / float(HASHTABLE_SIZE)
print '---'

worker_q = Queue.Queue()

for i in range(0, n_subsets):
    print '->', i
    start = divvy[i]
    end = divvy[i+1]
    worker_q.put((ht, i, start, end))

threads = []
for n in range(N_THREADS):
    t = threading.Thread(target=worker, args=(worker_q,))
    threads.append(t)
    t.start()

print 'started threads'

# wait for threads
for t in threads:
    t.join()
