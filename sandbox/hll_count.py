#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
## using a HyperLogLog counter and a bloom filter to count unique kmers,
## comparing results

import khmer
import sys
from screed.fasta import fasta_iter

from hyperloglog.hll import HyperLogLog


class new_hll_counter_py(object):

    def __init__(self, bits=8):
        # choose the precision by choosing how many estimators to track.
        self.bits = 8
        self.alpha = self._get_alpha(self.bits)
        self.num_bins = 1 << bits
        self.bit_bins = [1L << i for i in range(160 - self.bits + 1)]

        self.estimators = [0] * self.num_bins

    def _get_alpha(self, b):
        if not (4 <= b <= 16):
            raise ValueError("b=%d should be in range [4 : 16]" % b)

        if b == 4:
            return 0.673

        if b == 5:
            return 0.697

        if b == 6:
            return 0.709

        return 0.7213 / (1.0 + 1.079 / (1 << b))

    def estimate_cardinality(self):
        import math
        alpha = self.alpha
        bits = self.bits
        bins = self.estimators

        # harmonic mean
        E = alpha * float(len(bins) ** 2) / sum(math.pow(2.0, -x) for x in bins)

        if E <= 2.5 * bits:             # Small range correction
            V = bins.count(0)           # count number or registers equal to 0
            return bits * math.log(bins / float(V)) if V > 0 else E
        elif E <= float(1L << 160) / 30.0:
            return E
        else:
            return -(1L << 160) * math.log(1.0 - E / (1L << 160))

    # 'rho' function to calculate the bit pattern to watch (string of 0s)
    # here, 'rho' is the number of 0s to the left of the first 'accuracy' bits.
    def rho(self, w):
        import bisect
        r = len(self.bit_bins) - bisect.bisect_right(self.bit_bins, w)
        return r

    # to add a number into the counter:
    def add(self, num):
        from hashlib import sha1
        # take the hash of 'num'
        num = str(num)
        hash = long(sha1(num).hexdigest(), 16)

        # here, 'bin' is determined by the first 'bits' bits of hash
        bin = hash & ((1 << self.bits) - 1)

        # now count the number of 0s in the remaining bits
        remaining_bits = hash >> self.bits
        count = self.rho(remaining_bits)

        # take max of currently stored estimation & this one
        self.estimators[bin] = max(self.estimators[bin], count)


filename = sys.argv[1]
K = int(sys.argv[2])  # size of kmer
HT_SIZE = int(sys.argv[3])  # size of hashtable
N_HT = int(sys.argv[4])  # number of hashtables

ht = khmer.new_hashbits(K, HT_SIZE, N_HT)
hllcpp = khmer.new_hll_counter(8)
hllcpy = new_hll_counter_py(17)
hlllib = HyperLogLog(0.01)

n_unique = 0
for n, record in enumerate(fasta_iter(open(filename))):
    sequence = record['sequence']
    seq_len = len(sequence)
    for n in range(0, seq_len + 1 - K):
        kmer = sequence[n:n + K]
        if (not ht.get(kmer)):
            n_unique += 1
        ht.count(kmer)
        hllcpy.add(kmer)
        hllcpp.add(kmer)
        hlllib.add(kmer)

print 'unique:', n_unique
print 'bloom unique:', ht.n_unique_kmers()
print 'HLL py unique:', hllcpy.estimate_cardinality()
print 'HLL cpp unique:', hllcpp.estimate_cardinality()
print 'HLL lib unique:', len(hlllib)
