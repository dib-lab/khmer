## test the modified c++ code in hashbits.hh
# 1. n_occupied : faster counting of occupied bins in hashtables
# 2. n_unique_kmers : using bloom filter to count unique kmers ( use multiple hashtables)

import khmer
import sys
import screed
from screed.fasta import fasta_iter
import os
thisdir = os.path.dirname(__file__)
thisdir = os.path.abspath(thisdir)


def test_n_occupied_1():
    filename = os.path.join(thisdir, 'test-data/random-20-a.fa')
    #filename = 'test-data/random-20-a.fa'
    K = 20 # size of kmer
    HT_SIZE= 100000 # size of hashtable
    N_HT = 1 # number of hashtables

    ### test modified c++ n_occupied code
    ht1 = khmer.new_hashbits(K, HT_SIZE, N_HT)
    
    for n, record in enumerate(fasta_iter(open(filename))):
        ht1.consume(record['sequence'])
        
    assert ht1.n_occupied() == 3877  # this number is from original n_occupied code
    #print ht1.n_occupied()
def test_bloom_python_1():
        
    ### test python code to count unique kmers using bloom filter 
    filename = os.path.join(thisdir, 'test-data/random-20-a.fa')
    #filename = 'test-data/random-20-a.fa'
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
    #print n_unique
    print ht2.n_occupied()# 3882  ##### this should be 3882 since the process should have the equal effect as consume() method in test_bloom_c_1() function.
    assert  ht2.n_unique_kmers() == 3960 # this number equals to n_unique
    #print ht2.n_unique_kmers()

def test_bloom_c_1():
        
    ### test c++ code to count unique kmers using bloom filter 
    filename = os.path.join(thisdir, 'test-data/random-20-a.fa')
    #filename = 'test-data/random-20-a.fa'
    K = 20 # size of kmer
    HT_SIZE= 100000 # size of hashtable
    N_HT = 3 # number of hashtables
    
    ht3 = khmer.new_hashbits(K, HT_SIZE, N_HT)
    
    for n, record in enumerate(fasta_iter(open(filename))):
        ht3.consume(record['sequence'])
            
    assert ht3.n_occupied() == 3882
    assert  ht3.n_unique_kmers() == 3960
    #print ht3.n_unique_kmers()

def test_n_occupied_2(): # simple one
    K=4
    HT_SIZE = 10 # use 11
    N_HT = 1
    ht1 = khmer.new_hashbits(K, HT_SIZE, N_HT)
    ht1.count('AAAA') # 00 00 00 00 = 0
    assert ht1.n_occupied() == 1
    ht1.count('ACTG') # 00 10 01 11 = 
    assert ht1.n_occupied() == 2
    ht1.count('AACG') # 00 00 10 11 = 11  # collision 1
    #print ht1.n_occupied()
    assert ht1.n_occupied() == 2
    ht1.count('AGAC')   # 00  11 00 10 # collision 2
    assert ht1.n_occupied() == 2
#test_n_occupied_2()    

def test_bloom_c_2(): # simple one
    K=4
    HT_SIZE = 10 # use 11
    N_HT1 = 1   # hashtable size = 11
    N_HT2 = 2  # hashtable size = 11,13
    
    # use only 1 hashtable, no bloom filter
    ht1 = khmer.new_hashbits(K, HT_SIZE, N_HT1) 
    ht1.count('AAAA') # 00 00 00 00 = 0
    ht1.count('ACTG') # 00 10 01 11 =
    assert ht1.n_unique_kmers() == 2
    ht1.count('AACG') # 00 00 10 11 = 11  # collision  with 1st kmer
    assert ht1.n_unique_kmers() == 2
    ht1.count('AGAC')   # 00  11 00 10 # collision  with 2nd kmer
    assert ht1.n_unique_kmers() == 2

    # use two hashtables with 11,13
    ht2 = khmer.new_hashbits(K, HT_SIZE, N_HT2)
    ht2.count('AAAA') # 00 00 00 00 = 0

    ht2.count('ACTG') # 00 10 01 11 = 2*16 +4 +3 = 39 
    assert ht2.n_unique_kmers() == 2
    ht2.count('AACG') # 00 00 10 11 = 11  # collision with only 1st kmer
    assert ht2.n_unique_kmers() == 3
    ht2.count('AGAC')   # 00  11 00 10  3*16 +2 = 50 # collision with both 2nd and 3rd kmers
    assert ht2.n_unique_kmers() == 3
    
test_n_occupied_1()
test_n_occupied_2()
test_bloom_python_1()
test_bloom_c_1()
test_bloom_c_2()
