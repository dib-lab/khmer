import sys
import khmer

K=32
HASHTABLE_SIZE=int(116e9)
N_HT = 2

ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)

print 'loading mask from', sys.argv[1]
ht.consume_fasta(sys.argv[1])

print 'filtering', sys.argv[2]
ht.filter_if_present(sys.argv[2], sys.argv[2] + '.masked')
