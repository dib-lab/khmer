import khmer, sys
import os

K = 32
SUBSET_SIZE=int(1e4)
COUNTING_HT_SIZE=3e8                    # number of bytes
COUNTING_HT_N=4                         # number of counting hash tables
EXCURSION_DISTANCE=40
EXCURSION_KMER_THRESHOLD=200
EXCURSION_KMER_COUNT_THRESHOLD=5

##

filename = sys.argv[1]
basename = os.path.basename(filename)

counting = khmer.new_counting_hash(K, COUNTING_HT_SIZE, COUNTING_HT_N)

ht = khmer.new_hashbits(1, 1, 1)
print 'loading ht %s.ht' % basename
ht.load(basename + '.ht')

print 'loading tagset %s.tagset...' % basename
ht.load_tagset(basename + '.tagset')

divvy = ht.divide_tags_into_subsets(SUBSET_SIZE)
start, end = divvy[:2]
print 'doing pre-partitioning from', start, 'to', end

subset = ht.do_subset_partition(start, end)

ht.repartition_largest_partition(subset, counting,
                                 EXCURSION_DISTANCE,
                                 EXCURSION_KMER_THRESHOLD,
                                 EXCURSION_KMER_COUNT_THRESHOLD)

print 'saving stop tags'
ht.save_stop_tags(basename + '.stoptags')
