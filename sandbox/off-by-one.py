#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import sys
import screed
import khmer

x = [0] * 256
y = [0] * 256

ht = khmer.new_counting_hash(32, 3e8, 4)

for record in screed.open(sys.argv[1]):
    for i in range(0, len(record.sequence) - 31):
        kmer = record.sequence[i:i + 32]
        ht.count(kmer)

for record in screed.open(sys.argv[1]):
    for i in range(0, len(record.sequence) - 31):
        kmer = record.sequence[i:i + 32]
        n = ht.get(kmer)
        m = ht.max_hamming1_count(kmer)

        x[n] += 1
        y[m] += 1

for i, (n, m) in enumerate(zip(x, y)):
    print "%d,%d,%d" % (i, n, m,)
