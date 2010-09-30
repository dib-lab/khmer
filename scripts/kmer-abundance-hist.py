#! /usr/bin/env python
import sys, khmer

K = 32

###

output = sys.argv[1]
fa_files = sys.argv[2:]

HT_SIZE = int(1e9)
HT_SIZE = khmer.get_n_primes_above_x(1, HT_SIZE)[0]
print HT_SIZE

ht = khmer.new_hashtable(K, HT_SIZE)

for filename in fa_files:
    ht.consume_fasta(filename)

print 'preparing hist...'
z = ht.abundance_distribution()
fp = open(output, 'w')
for n, i in enumerate(z[1:]):
    print >>fp, n + 1, i
