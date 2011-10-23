import khmer, sys, os
import gc
import glob

already_part = sys.argv[1]
new_to_part = sys.argv[2]
basename = os.path.basename(new_to_part)
pmap_filename = sys.argv[3]

if not os.path.exists(already_part):
    print '%s doesn\'t exist! dying.' % already_part
    sys.exit(0)

# create a fake-ish ht; K matters, but not hashtable size.
ht = khmer.load_hashbits(already_part + '.ht')
ht.load_tagset(already_part + '.tagset')
ht.merge_subset_from_disk(pmap_filename)

# find singletons
n_singletons = ht.find_unpart(new_to_part)
print 'found:', n_singletons

n_partitions = ht.output_partitions(new_to_part, basename + '.unpart')

###

(n_partitions, n_singletons) = ht.count_partitions()

print 'output partitions:', n_partitions
print 'pmap partitions:', n_partitions
print 'singletons:', n_singletons
