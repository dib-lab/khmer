#!/usr/bin/env python

# A demonstration of using khmer for exact k-mer counting. The memory required
# is 4^k, which limits this to small values of k.

import khmer

# Note:
#    - The forward and reverse complements will be collapsed since in this case
#      k is even.
#    - There are 4^k possible sequences of length k.
#    - If the table size provided to the countgraph is not a prime number, it
#      will select the next lowest prime number. So here we are requesting a
#      table size of *slightly more* than 4^k rather than *slightly less* so we
#      can avoid any false positives.
ksize = 6
nkmers = 4**ksize
tablesize = nkmers + 10

# Initialize countgraph
cg = khmer.Countgraph(ksize, tablesize, 1)
print('Created a countgraph with', cg.hashsizes(), 'buckets')

# Increment the count of some k-mers
cg.count('ATGGCA')
cg.count('ATGGCA')
cg.count('ACATGG')
cg.count('AAAAAA')
cg.count('TTTTTT')  # this will be counted towards AAAAAA

# Show all >0 k-mer abundances from the table
for i in range(nkmers):
    if cg.get(i):
        print(cg.reverse_hash(i), cg.get(i))


# Note: The reverse_hash function is only available for Countgraph and
# Nodegraph, not Counttable and Nodetable.
