#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# make overlapping & connectable reads from a random DNA sequence, for testing
# purposes.
import random

K = 31
SIZE = 100000

dna = ["A"] * SIZE + ["C"] * SIZE + ["T"] * SIZE + ["G"] * SIZE
random.shuffle(dna)

n = 0
dna = "".join(dna)

x = []
for i in range(40, len(dna), 40):
    subseq = dna[i - (K - 1):i + 40]
    x.append((n, subseq))
    n += 1

random.shuffle(x)
for n, subseq in x:
    print '>%d\n%s' % (n, subseq)
