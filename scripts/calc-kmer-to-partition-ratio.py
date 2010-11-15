#! /usr/bin/env python
import sys
import khmer
from screed.fasta import fasta_iter

K = 32
HASHTABLE_SIZE=int(8e9)
N_HASHTABLES=4

total_kmers = 0

ht = khmer.new_hashbits(32, HASHTABLE_SIZE, N_HASHTABLES)
pidset = set()
for n, record in enumerate(fasta_iter(open(sys.argv[1]), parse_description=False)):
    pid = record['name'].rsplit('\t', 1)[1]
    pidset.add(pid)
    
    ht.consume(record['sequence'])
    total_kmers += len(record['sequence']) - K + 1

print 'n partitions:', len(pidset)
print 'unique kmers:', ht.n_unique_kmers()
print 'total kmers:', total_kmers
