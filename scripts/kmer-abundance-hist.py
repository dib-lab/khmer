#! /usr/bin/env python
import sys, khmer

K = int(sys.argv[1])
output = sys.argv[2]
fa_files = sys.argv[3:]

HT_SIZE = min(4**K, 4**17+1)

ht = khmer.new_hashtable(K, HT_SIZE)

for filename in fa_files:
    ht.consume_fasta(filename)

z = ht.abundance_distribution()
fp = open(output, 'w')
for n, i in enumerate(z[1:]):
    print >>fp, n + 1, i
