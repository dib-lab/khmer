#! /usr/bin/env python
import khmer, sys, screed

K=32
HASHTABLE_SIZE=int(4**15)+1

infile = sys.argv[1]
outfile = sys.argv[2]
min_partition_size = int(sys.argv[3])

ht = khmer.new_hashtable(K, HASHTABLE_SIZE)

for n, record in enumerate(screed.fasta.fasta_iter(open(infile))):
    if n % 10000 == 0:
        print '...', n
    seq = record['sequence']
    try:
        ppi = ht.consume_and_find_tags(seq)
        if ppi:
            ht.assign_partition_id(ppi)
    except ValueError:
        pass

n_kept = ht.output_partitions(infile, outfile)
print n_kept, 'partitions kept'
