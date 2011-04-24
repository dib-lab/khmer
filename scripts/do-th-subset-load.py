import khmer, sys, os
import gc
import glob

K = 32

filename=sys.argv[1]
basename = os.path.basename(filename)
subset_filenames=sys.argv[2:]

if not os.path.exists(filename):
    print '%s doesn\'t exist! dying.' % filename
    sys.exit(0)

print '---'
print 'outputting partitioned file to', basename + '.part'
print 'merged pmap will be in', basename + '.pmap.merged'
print '---'

# create a fake-ish ht; K matters, but not hashtable size.
ht = khmer.new_hashbits(K, 1, 1)
ht.load(basename + '.ht')
ht.load_tagset(basename + '.tagset')
 
# load & merge
counting = khmer.new_counting_hash(K, 3e8, 4)
for subset_file in subset_filenames:
    print '<-', subset_file
    subset = ht.load_subset_partitionmap(subset_file)

    last_size = 0
    size = 1001
    while size > 1000 and last_size != size:
        last_size = size
        size = ht.repartition_largest_partition(subset, counting, 40, 200, 3)
        print 'repart40:', subset_file

    if 1:
        last_size = 0
        size = 1001
        while size > 1000 and last_size != size:
            last_size = size
            size = ht.repartition_largest_partition(subset, counting, 1000, 5000, 5)
            print 'repart100:', subset_file

    #os.rename(subset_file, subset_file + '.bak')
    #ht.save_subset_partitionmap(subset, subset_file + '.new')

    print 'saving stoptags binary'
    ht.save_stop_tags(basename + '.stoptags')
