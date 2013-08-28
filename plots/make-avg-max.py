#! /usr/bin/env python
#
# This script is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import sys

cutoff = int(sys.argv[2])

group_l = {}

for line in open(sys.argv[1]):
    _, group_n, size = line.split()
    group_n = int(group_n)
    size = int(size)
    if size >= cutoff:
        x = group_l.get(group_n, [])
        x.append(size)
        group_l[group_n] = x

max_g = {}
avg_g = {}

for group_n in group_l:
    sizes = group_l[group_n]
    max_g[group_n] = max(sizes)
    avg_g[group_n] = sum(sizes) / float(len(sizes))

for group_n in sorted(group_l):
    print group_n, len(group_l[group_n]), avg_g[group_n], max_g[group_n]
