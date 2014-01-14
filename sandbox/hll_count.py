#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
## using a HyperLogLog counter and a bloom filter to count unique kmers,
## comparing results

from collections import Counter

import khmer
import sys
from screed.fasta import fasta_iter

from hyperloglog.hll import HyperLogLog


filename = sys.argv[1]
K = int(sys.argv[2])  # size of kmer
HT_SIZE = int(sys.argv[3])  # size of hashtable
N_HT = int(sys.argv[4])  # number of hashtables

ht = khmer.new_hashbits(K, HT_SIZE, N_HT)
hllcpp = khmer.new_hll_counter(0.01)
hlllib = HyperLogLog(0.01)
counter = Counter()

n_unique = 0
for n, record in enumerate(fasta_iter(open(filename))):
    sequence = record['sequence']
    seq_len = len(sequence)
    for n in range(0, seq_len + 1 - K):
        kmer = sequence[n:n + K]
#        if (not ht.get(kmer)):
#            n_unique += 1
#        ht.count(kmer)
        counter.update([kmer])
        hllcpp.add(kmer)
        hlllib.add(kmer)

print 'unique:', n_unique
print 'bloom unique:', ht.n_unique_kmers()
print 'HLL cpp unique:', hllcpp.estimate_cardinality()
print 'HLL lib unique:', len(hlllib)
print 'Python stdlib counter:', len(counter)
