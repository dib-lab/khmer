#! /usr/bin/env python
import sys
import khmer
import gc

MIN_PARTITION_SIZE = 3
K = 32

tagset_filename = sys.argv[1]
subset_filenames = sys.argv[2:]

# 
ht = khmer.new_hashtable(32, 1)
tagmap = ht.load_tagmap(tagset_filename)

for filename in subset_filenames:
    print 'maxifying:', filename
    subset = ht.load_subset_partitionmap(filename)
    ht.subset_maxify_partition_size(subset, tagmap)
    del subset
    gc.collect()

print 'discarding'
ht.discard_tags(tagmap, MIN_PARTITION_SIZE)

for filename in subset_filenames:
    print 'filtering:', filename
    subset = ht.load_subset_partitionmap(filename)
    ht.subset_filter_against_tags(subset, tagmap)
    
    new_filename = 'filtered_' + filename[:-5]
    ht.save_subset_partitionmap(subset, new_filename + '.pmap')

    del subset
    gc.collect()
