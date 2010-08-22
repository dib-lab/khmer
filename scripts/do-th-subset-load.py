import khmer, sys
import threading
ht = khmer.new_hashtable(32, 1)

filename=sys.argv[1]
SUBSET_SIZE = 100000

def load(filename, ht, start, stop):
    outfile = filename + '.subset.%d-%d' % (start, stop)
    subset = ht.load_subset_partitionmap(outfile + '.pmap', outfile + '.surr')
    ht.merge_subset(subset)

(total_reads, total_kmers) = ht.consume_fasta_and_tag(filename)
n_subsets = total_reads / SUBSET_SIZE + 1

for i in range(0, n_subsets):
    print '<-', i
    start = i*SUBSET_SIZE
    stop = (i + 1)*SUBSET_SIZE

    load(filename, ht, start, stop)

n_partitions = ht.output_partitions(filename, filename + '.part')
print ht.count_partitions()
