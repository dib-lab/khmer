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
            (ht, filename, start, stop) = q.get(False)
        except Queue.Empty:
            print 'exiting'
            return

        outfile = filename + '.subset.%d-%d' % (start, stop)
        if os.path.exists(outfile + '.pmap'):
            print 'SKIPPING', filename, start, stop
            continue
        
        print 'starting:', filename, start, stop
        subset = ht.do_subset_partition(filename, start, stop)
        print 'saving:', filename, start, stop
        
        outfile = filename + '.subset.%d-%d' % (start, stop)
        ht.save_subset_partitionmap(subset,
                                    outfile + '.pmap',
                                    outfile + '.surr')
        del subset
        gc.collect()

(total_reads, total_kmers) = ht.consume_fasta_and_tag(FILENAME)
n_subsets = total_reads / SUBSET_SIZE + 1

print '---'
print 'hashtable occupancy:', ht.n_occupied() / float(HASHTABLE_SIZE)
print '---'

worker_q = Queue.Queue()

for i in range(0, n_subsets):
    print '->', i
    start = i*SUBSET_SIZE
    end = (i + 1)*SUBSET_SIZE
    worker_q.put((ht, FILENAME, start, end))

threads = []
for n in range(N_THREADS):
    t = threading.Thread(target=worker, args=(worker_q,))
    threads.append(t)
    t.start()

print 'started threads'

# wait for threads
for t in threads:
    t.join()
