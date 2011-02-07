import khmer, sys
import threading
import Queue
import gc
import os.path

K=32
HASHTABLE_SIZE=int(8e9)
N_HT=4

SUBSET_SIZE = int(2e5)
N_THREADS = 8

ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)

###

save_ht = False
load_ht = False
save_merged_pmap = True

assert not (save_ht and load_ht)         # incompatible

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
        subset = ht.do_subset_partition(start, stop)

        print 'saving:', basename, n
        ht.save_subset_partitionmap(subset, outfile)
        del subset
        gc.collect()

def main(filename):
    global ht
    
    basename = os.path.basename(filename)

    print 'K', K
    print 'HASHTABLE SIZE %g' % HASHTABLE_SIZE
    print 'SUBSET SIZE', SUBSET_SIZE
    print 'N THREADS', N_THREADS
    print '--'

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
        ht.save(basename + '.ht')
        print 'loading tagset %s.tagset...' % basename
        ht.save_tagset(basename + '.tagset')

    # divide the tags up into subsets
    divvy = ht.divide_tags_into_subsets(SUBSET_SIZE)
    n_subsets = len(divvy)
    divvy.append(0)

    # build a queue of tasks:
    worker_q = Queue.Queue()

    for i in range(0, n_subsets):
        start = divvy[i]
        end = divvy[i+1]
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
    
    # create a new, empty ht object for merging; K matters, but not
    # hashtable size.
    ht = khmer.new_hashbits(32, 1, 1)

    # load & merge all pmap files
    for i in range(0, n_subsets):
        pmap_file = basename + '.subset.%d.pmap' % (i,)
        print 'loading', pmap_file
        ht.merge_subset_from_disk(pmap_file)

    # save merged partitionmap
    if save_merged_pmap:
        print 'saving merged pmap to %s.pmap.merged' % basename 
        ht.save_partitionmap(basename + '.pmap.merged')

    # output partitions!
    n_partitions = ht.output_partitions(filename, basename + '.part')
    (n_partitions, n_singletons) = ht.count_partitions()
    print 'output partitions:', n_partitions
    print 'pmap partitions:', n_partitions
    print 'singletons:', n_singletons
    

if __name__ == '__main__':
    main(sys.argv[1])
