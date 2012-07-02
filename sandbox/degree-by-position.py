#! /usr/bin/env python
import khmer, sys
from screed.fasta import fasta_iter

K = 32
HASHTABLE_SIZE=int(8e9)
N_HT = 4

outfp = open(sys.argv[2], 'w')

ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)
ht.consume_fasta(sys.argv[1])

hist = [0] * 200
histcount = [0] * 200
for n, record in enumerate(fasta_iter(open(sys.argv[1]))):
    if n % 10000 == 0:
        print '...', n

    seq = record['sequence']
    for pos in range(0, len(seq) - K + 1):
        kmer = seq[pos:pos+K]
        count = ht.kmer_degree(kmer)

        hist[pos] += count
        histcount[pos] += 1

for i in range(len(hist)):
    total = hist[i]
    count = histcount[i]
    if not count:
        continue
    
    print >>outfp, i, total, count, total / float(count)
