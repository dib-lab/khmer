#! /usr/bin/env python
import sys, khmer

K = 32
N_HT=4
HT_SIZE = int(5e6)

###

hashfile = sys.argv[1]
filename = sys.argv[2]
output = sys.argv[3]

ht = khmer.load_counting_hash(hashfile)
tracking = khmer.new_hashbits(K, HT_SIZE, N_HT)

print 'preparing hist...'
z = ht.abundance_distribution(filename, tracking)

fp = open(output, 'w')
for n, i in enumerate(z[1:]):
    print >>fp, n + 1, i
