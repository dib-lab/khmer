import khmer, sys
import threading

K = 32
HASHTABLE_SIZE=22906493
N_THREADS=4                             # @CTB

ht = khmer.new_hashbits(K, HASHTABLE_SIZE)

###

def threaded_calc(ht, start, stop, results):
    print 'starting', start, stop
    x = ht.do_subset_partition(start, stop)
    print 'done!', start, stop
    results.append(x)

def main(filename):
    print 'K', K
    print 'HASHTABLE SIZE %g' % HASHTABLE_SIZE
    print 'N THREADS', N_THREADS
    print '--'

    (total_reads, total_kmers) = ht.consume_fasta_and_tag(filename)
    print total_reads

    subset_size = total_reads / N_THREADS + total_reads % N_THREADS

    results = []

    #ht.save(filename + '.ht')                # @CTB
    #ht.save_tagset(filename + '.tagset')     # @CTB
    #ht.load(filename + '.ht')
    #ht.load_tagset(filename + '.tagset')

    # calculate the hashtable occupancy
    print '---'
    #print 'hashtable occupancy:', ht.n_occupied() / float(HASHTABLE_SIZE)
    print '---'

    divvy = ht.divide_tags_into_subsets(subset_size)
    n_subsets = len(divvy)
    divvy.append(0)
        
    threads = []

    # start things.
    for i in range(n_subsets):
        start = divvy[i]
        end = divvy[i+1]
        t = threading.Thread(target=threaded_calc,
                             args=(ht, start, end, results))
        threads.append(t)
        t.start()

    print 'started:', N_THREADS - 1

    # wait for them all to end.
    for i, t in enumerate(threads):
        t.join()
        print 'done: ', i

    # merge
    for i, x in enumerate(results):
        print 'merging %d' % (i,)
        ht.merge_subset(x)
    
    n_partitions = ht.output_partitions(filename, filename + '.part')
    print n_partitions, 'partitions kept'
    print 'n surrendered:', ht.count_partitions()[2]

if __name__ == '__main__':
    main(sys.argv[1])
