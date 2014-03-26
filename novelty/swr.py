#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import khmer
import random
import sys

K = 12                                  # size of K
N = 25000                               # 1/4 the size of the genome
P_ERROR = .01                           # per-base probability of error

###

# construct a random genome
genome = "A"*N + "C"*N + "G"*N + "T"*N
genome = list(genome)
random.shuffle(genome)
genome = "".join(genome)

###

# count the number of unique k-mers
kt = khmer.new_ktable(K)
kt.consume(genome)

total = 0
for i in range(0, 4**K):
    if kt.get(i):
        total += 1

print >> sys.stderr, "%d unique k-mers in genome" % total

###

# go through, sample with replacement and mutation, and calculate
# number of novel k-mers picked as a function of sampling.

kt = khmer.new_ktable(K)
n = 0
for i in range(0, 20*4*N):
    # pick random k-mer
    pos = random.randint(0, 4*N - K)
    subseq = genome[pos:pos+K]

    # should we mutate?
    if random.uniform(0, 1) <= P_ERROR:
        z = random.choice(range(K))
        subseq = list(subseq)
        new_bp = random.choice('ACGT')
        while subseq[z] == new_bp:
            new_bp = random.choice('ACGT')
            
        subseq[z] = new_bp
        subseq = "".join(subseq)

    # count!
    kt.count(subseq)

    # is it novel?
    if kt.get(subseq) == 1:
        n += 1

    # progress report and output
    if i % 1000 == 0:
        if i % 100000 == 0:
            print>>sys.stderr, '...', i
        print i, n
