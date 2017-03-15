#!/usr/bin/env python

# A demonstration of using khmer for exact k-mer counting. The memory required
# is 4^k, which limits this to small values of k.

from __future__ import print_function
import khmer

# Note: forward and reverse complement will be collapsed since k is even.
ksize = 6
nkmers = 4**ksize
cg = khmer.Countgraph(ksize, nkmers, 1)

cg.count('ATGGCA')
cg.count('ATGGCA')
cg.count('ACATGG')
cg.count('AAAAAA')
cg.count('TTTTTT') # this will be counted towards AAAAAA

for i in range(0, nkmers):
    if cg.get(i):
        print(cg.reverse_hash(i), cg.get(i))
