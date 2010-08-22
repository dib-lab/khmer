import khmer, sys
import threading
ht = khmer.new_hashtable(32, 4**15+1)

FILENAME=sys.argv[1]
N_THREADS=8

(total_reads, total_kmers) = ht.consume_fasta_and_tag(FILENAME)
print total_reads

subset_size = total_reads / N_THREADS

results = []
def threaded_calc(ht, filename, start, stop):
    x = ht.do_subset_partition(filename, start, stop)
    print 'done!'
    results.append(x)
        
threads = []

# start things.
for i in range(N_THREADS - 1):
    start = subset_size * i
    end = start + subset_size
    t = threading.Thread(target=threaded_calc, args=(ht, FILENAME, start, end))
    threads.append(t)
    
    print 'starting', i
    t.start()

start = subset_size * (N_THREADS - 1)
end = total_reads
t = threading.Thread(target=threaded_calc, args=(ht, FILENAME, start, end))
threads.append(t)
t.start()
print 'starting', N_THREADS - 1

# wait for them all to end.
for i, t in enumerate(threads):
    t.join()
    print 'done: ', i

# merge
for i, x in enumerate(results):
    print 'merging %d' % (i,)
    ht.merge_subset(x)
    
n_partitions = ht.output_partitions(FILENAME, FILENAME + '.part')
print n_partitions, 'partitions kept'
