#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
## using a HyperLogLog counter and a bloom filter to count unique kmers,
## comparing results

from collections import Counter
import string

import khmer
import sys
from screed.fasta import fasta_iter

from hyperloglog.hll import HyperLogLog


filename = sys.argv[1]
K = int(sys.argv[2])  # size of kmer

ERROR_RATE = .01
TT = string.maketrans('ACGT', 'TGCA')

hllcpp = khmer.new_hll_counter(ERROR_RATE)
hlllib = HyperLogLog(ERROR_RATE)
counter = Counter()
counter_norc = Counter()

for n, record in enumerate(fasta_iter(open(filename))):
    sequence = record['sequence']
    seq_len = len(sequence)
    for n in range(0, seq_len + 1 - K):
        kmer = sequence[n:n + K]
        rc = kmer[::-1].translate(TT)

        hllcpp.add(kmer)
        hlllib.add(kmer)
        counter_norc.update([kmer])

        if rc in counter:
            kmer = rc
        counter.update([kmer])
    #hllcpp.consume_string(sequence, K)

cpp_estimate = hllcpp.estimate_cardinality()
py_estimate = len(hlllib)
real_count = len(counter)

print 'Unique:', real_count
print 'Unique no rc:', len(counter_norc)

print 'HLL cpp unique:', cpp_estimate,
print ', error:', round(float(abs(cpp_estimate - real_count)) / (real_count or 1), 3)

print 'HLL lib unique:', py_estimate,
print ', error:', round(float(abs(py_estimate - real_count)) / (real_count or 1), 3)
print
