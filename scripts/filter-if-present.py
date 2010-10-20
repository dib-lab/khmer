import sys
import khmer

K=32
HASHTABLE_SIZE=int(1e9)
N_HT = 2

ht = khmer.new_hashbits(32, 8e9, 2)

print 'loading mask from', sys.argv[1]
ht.consume_fasta(sys.argv[1])

print 'filtering', sys.argv[2]
ht.filter_if_present(sys.argv[2], sys.argv[2] + '.masked')
