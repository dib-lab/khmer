import sys, screed.fasta

K = 32
HASHTABLE_SIZE=int(8e9)
N_HT = 4
THRESHOLD=100

import khmer
ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)

ht.consume_fasta(sys.argv[1])

for record in screed.fasta.fasta_iter(open(sys.argv[1])):
    seq = record['sequence']
    print ht.count_kmers_within_radius(seq[:K], THRESHOLD)
