import khmer, sys
import gc
import glob

K = 32

filename=sys.argv[1]
subset_filenames=sys.argv[2:]

def load(filename, ht):
    pmap_filename = filename
    print 'loading', filename
    subset = ht.load_subset_partitionmap(pmap_filename)
    print 'merging', filename
    ht.merge_subset(subset)

# create a fake-ish ht; K matters, but not hashtable size.
ht = khmer.new_hashtable(32, 1)

# detect all of the relevant partitionmap files
#subset_filenames = glob.glob(filename + '.subset.*.pmap')

# load & merge
for subset_file in subset_filenames:
    print '<-', subset_file
    load(subset_file, ht)

# partition!
n_partitions = ht.output_partitions(filename, filename + '.part')
print n_partitions
print ht.count_partitions()
