import khmer, sys
import threading
ht = khmer.new_hashtable(32, 4**15+11)

FILENAME=sys.argv[1]
N_THREADS=8

(total_reads, total_kmers) = ht.consume_fasta_and_tag(FILENAME)
print total_reads

subset_size = total_reads / N_THREADS + total_reads % N_THREADS

results = []
def threaded_calc(ht, start, stop):
    x = ht.do_subset_partition(start, stop)
    print 'done!'
    results.append(x)

divvy = ht.divide_tags_into_subsets(subset_size)
n_subsets = len(divvy)
divvy.append(0)
        
threads = []
# start things.
for i in range(n_subsets):
    start = divvy[i]
    end = divvy[i+1]
    t = threading.Thread(target=threaded_calc, args=(ht, start, end))
    threads.append(t)
    
    print 'starting', i
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
print 'n surrendered:', ht.count_partitions()[2]

