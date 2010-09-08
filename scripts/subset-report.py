import khmer, sys
import gc
import glob

K = 32

subset_filenames=sys.argv[1:]

ht = khmer.new_hashtable(K, 1)
for filename in subset_filenames:
    print 'partition map:', filename
    subset = ht.load_subset_partitionmap(filename)
    n_part, n_orphan, n_surr = ht.subset_count_partitions(subset)
    print 'num partitions:', n_part
    print 'num orphans:', n_orphan
    print 'surrendered:', n_surr
    print '--'
