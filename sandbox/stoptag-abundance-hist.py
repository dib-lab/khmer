#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import sys
import khmer
import os
from pylab import *

hashfile = sys.argv[1]
filename = sys.argv[2]
figure = sys.argv[3]

ht = khmer.load_counting_hash(hashfile)

outabund = open(os.path.basename(filename) + '.counts', 'w')

counts = []
d = {}
for sequence in open(sys.argv[2]):
    sequence = sequence.strip()

    count = ht.get(sequence)
    counts.append(count)
    d[count] = d.get(count, 0) + 1

    if count > 1000:
        print >>outabund, sequence, count

outfp = open(figure + '.countshist', 'w')
sofar = 0
sofar_cumu = 0
for k in sorted(d.keys()):
    sofar += d[k]
    sofar_cumu += k * d[k]
    print >>outfp, k, d[k], sofar, sofar_cumu

hist(counts, normed=True, cumulative=True, bins=100, range=(1, 1000))
savefig(figure)
