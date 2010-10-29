## using bloom filter to count unique kmers

import khmer
import sys
import screed
from screed.fasta import fasta_iter

filename = sys.argv[1]
K = int(sys.argv[2]) # size of kmer
HT_SIZE= int(sys.argv[3])# size of hashtable
N_HT = int(sys.argv[4]) # number of hashtables

ht = khmer.new_hashbits(K, HT_SIZE, N_HT)

n_unique = 0
for n, record in enumerate(fasta_iter(open(filename))):
    sequence = record['sequence']
    seq_len = len(sequence)
    for n in range(0,seq_len+1-K):
        kmer = sequence[n:n+K]
        if (not ht.get(kmer)):
            n_unique+=1
        ht.count(kmer)
        
print n_unique
print ht.n_occupied()
print ht.n_unique_kmers()
