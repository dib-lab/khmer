import khmer, sys
import gc
import glob

K = 32

filename=sys.argv[1]
subset_filenames=sys.argv[2:]

def load(filename, ht):
    print 'loading', filename
    subset = ht.merge_subset_from_disk(filename)

# create a fake-ish ht; K matters, but not hashtable size.
ht = khmer.new_hashbits(32, 1, 1)

# load & merge
for subset_file in subset_filenames:
    print '<-', subset_file
    load(subset_file, ht)

# partition!
n_partitions = ht.output_partitions(filename, filename + '.part')
print n_partitions
print ht.count_partitions()
