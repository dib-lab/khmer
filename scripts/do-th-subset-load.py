import khmer, sys, os
import gc
import glob

K = 32

filename=sys.argv[1]
subset_filenames=sys.argv[2:]

if not os.path.exists(filename):
    print '%s doesn\'t exist! dying.' % filename
    sys.exit(0)

# create a fake-ish ht; K matters, but not hashtable size.
ht = khmer.new_hashbits(32, 1, 1)
 
# load & merge
for subset_file in subset_filenames:
    print '<-', subset_file
    ht.merge_subset_from_disk(subset_file)

# save merged partitionmap
ht.save_partitionmap(filename + '.pmap.merged')

# partition!
n_partitions = ht.output_partitions(filename, filename + '.part')
print n_partitions
print ht.count_partitions()
