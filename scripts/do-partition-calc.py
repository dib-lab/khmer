#! /usr/bin/env python
import khmer, sys

K=32
HASHTABLE_SIZE=int(4**15)+1

infile = sys.argv[1]
outfile = sys.argv[2]
min_partition_size = int(sys.argv[3])

ht = khmer.new_hashtable(K, HASHTABLE_SIZE)

#n_partitions = ht.do_exact_partition(infile)
n_partitions = ht.do_truncated_partition(infile, outfile)
print n_partitions, 'partitions kept'
