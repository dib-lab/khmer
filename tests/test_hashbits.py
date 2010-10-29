## using bloom filter to count unique kmers

import khmer
import sys
import screed
from screed.fasta import fasta_iter

filename = 'test-data/random-20-a.fa'
K = 20 # size of kmer
HT_SIZE= 100000 # size of hashtable
N_HT = 1 # number of hashtables

### test modified c++ n_occupied code
ht1 = khmer.new_hashbits(K, HT_SIZE, N_HT)

for n, record in enumerate(fasta_iter(open(filename))):
    ht1.consume(record['sequence'])
    
assert ht1.n_occupied() == 3877  # this number is from original n_occupied code


### test python code to count unique kmers using bloom filter 
filename = 'test-data/random-20-a.fa'
K = 20 # size of kmer
HT_SIZE= 100000 # size of hashtable
N_HT = 3 # number of hashtables

ht2 = khmer.new_hashbits(K, HT_SIZE, N_HT)

n_unique = 0
for n, record in enumerate(fasta_iter(open(filename))):
    sequence = record['sequence']
    seq_len = len(sequence)
    for n in range(0,seq_len+1-K):
        kmer = sequence[n:n+K]
        if (not ht2.get(kmer)):
            n_unique+=1
        ht2.count(kmer)

assert n_unique == 3960
assert ht2.n_occupied() == 3882
assert  ht2.n_unique_kmers() == 3960 # this number equals to n_unique


### test c++ code to count unique kmers using bloom filter 
filename = 'test-data/random-20-a.fa'
K = 20 # size of kmer
HT_SIZE= 100000 # size of hashtable
N_HT = 3 # number of hashtables

ht3 = khmer.new_hashbits(K, HT_SIZE, N_HT)

for n, record in enumerate(fasta_iter(open(filename))):
    ht3.consume(record['sequence'])
        
assert ht3.n_occupied() == 3882
assert  ht3.n_unique_kmers() == 3960
