#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
import sys

cutoff = int(sys.argv[2])

sizes = []
for line in open(sys.argv[1]):
    _, _, size = line.split()
    size = int(size)
    if size >= cutoff:
        sizes.append(size)

sizes.sort()
sizes.reverse()

fp = open(sys.argv[1] + '.rank', 'w')
total = 0
d = {}
for rank, size in enumerate(sizes):
    total += size
    d[size] = total
    print >>fp, rank, size, total

totalsize = sum(sizes)

fp = open(sys.argv[1] + '.size', 'w')
total = 0
keys = d.keys()
keys.sort()
for k in keys:
    print >>fp, k, d[k], totalsize - d[k]

