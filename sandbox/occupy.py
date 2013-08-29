#! /usr/bin/env
import sys
import khmer

filename = sys.argv[1]
k = int(sys.argv[2])

ht = khmer.new_hashtable(k, 4 ** k)
ht.consume_fasta(filename)

N = ht.n_occupied()

for i in range(0, 4):
    print 'q', i, ht.n_occupied(i * 4 ** (k - 1), (i + 1) * 4 ** (k - 1))
