import khmer
import sys
import screed


input_fasta = sys.argv[3]
K = sys.argv[1]
x = sys.argv[2]


ht = khmer.new_hashbits(K, x, 4)


