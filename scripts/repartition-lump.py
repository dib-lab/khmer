import khmer, sys, os
import gc
import glob

K = 32

# counting hash parameters.
COUNTING_HT_SIZE=3e8                    # number of bytes
COUNTING_HT_N=4                         # number of counting hash tables

# Lump removal parameters.  Probably shouldn't be changed, but who knows?
#
# explanation:
#
# We will walk EXCURSION_DISTANCE out from each tag; if we find more than
# EXCURSION_KMER_THRESHOLD kmers within that range, this will be a "big"
# excursion and we will track all k-mers visited.  If we find that any
# k-mer has been visited more than EXCURSION_KMER_COUNT_THRESHOLD times,
# we will mark it as BAD and make it a stop tag for traversal.

EXCURSION_DISTANCE=40
EXCURSION_KMER_THRESHOLD=200
EXCURSION_KMER_COUNT_THRESHOLD=5

###

filename=sys.argv[1]
basename = os.path.basename(filename)
subset_filenames=sys.argv[2:]

if not os.path.exists(basename + '.ht'):
    print '%s.ht doesn\'t exist! dying.' % basename
    sys.exit(0)
    
if not os.path.exists(basename + '.tagset'):
    print '%s.tagset doesn\'t exist! dying.' % basename
    sys.exit(0)

print '---'
print 'stoptags will be in', basename + '.stoptags'
print '---'

ht = khmer.new_hashbits(K, 1, 1)
ht.load(basename + '.ht')
ht.load_tagset(basename + '.tagset')

# create counting hash
counting = khmer.new_counting_hash(K, COUNTING_HT_SIZE, COUNTING_HT_N)
 
# load & merge
for subset_file in subset_filenames:
    print '<-', subset_file
    subset = ht.load_subset_partitionmap(subset_file)

    print '** repartitioning subset... %s' % subset_file
    ht.repartition_largest_partition(subset, counting,
                                     EXCURSION_DISTANCE,
                                     EXCURSION_KMER_THRESHOLD,
                                     EXCURSION_KMER_COUNT_THRESHOLD)
    
    print '** merging subset... %s' % subset_file
    ht.merge_subset(subset)
    
    print '** repartitioning, round 2... %s' % subset_file
    size = ht.repartition_largest_partition(None, counting,
                                            EXCURSION_DISTANCE,
                                            EXCURSION_KMER_THRESHOLD,
                                            EXCURSION_KMER_COUNT_THRESHOLD)

    print '** repartitioned size:', size

    print 'saving stoptags binary'
    ht.save_stop_tags(basename + '.stoptags')
