#! /usr/bin/env python
#
# This script is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import sys
from screed.fasta import fasta_iter

MAX_SIZE = 5000

distfp = open(sys.argv[2], 'w')
count = {}
filename = sys.argv[1]

for n, record in enumerate(
        fasta_iter(open(sys.argv[1]), parse_description=False)):
    if n % 10000 == 0:
        print '...', n

    partition = int(record['name'].rsplit('\t', 1)[1])
    count[partition] = count.get(partition, 0) + 1

# develop histogram of partition sizes
dist = {}
for n, record in enumerate(
        fasta_iter(open(sys.argv[1]), parse_description=False)):
    if n % 10000 == 0:
        print '...x2', n

    partition = int(record['name'].rsplit('\t', 1)[1])
    if partition not in count:
        continue

    c = count[partition]
    if partition == 0:
        c = 0

    dist[c] = dist.get(c, 0) + 1

#    if c >= threshold:
#        outfp.write('>%s\n%s\n' % (record['name'], record['sequence']))

# output histogram
total = 0
for c, n in sorted(dist.items()):
    if c:
        n /= c
    total += n
    distfp.write('%d %d %d\n' % (c, n, total))

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
