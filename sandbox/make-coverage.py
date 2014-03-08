#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
import screed

import sys

dbfile = sys.argv[1]
mapfile = sys.argv[2]

lengths = {}
for n, record in enumerate(screed.open(dbfile)):
    if n % 10000 == 0:
        print '...', n
    lengths[record.name] = len(record.sequence)

sums = {}
for n, line in enumerate(open(mapfile)):
    if n % 10000 == 0:
        print '... 2x', n
    x = line.split('\t')
    name = x[2]
    readlen = len(x[4])
    sums[name] = sums.get(name, 0) + 1

mapped_reads = n

rpkms = {}
for k in sums:
    rpkms[k] = sums[k] * (1000. / float(lengths[k])) * float(mapped_reads)/1e6

outfp = open(dbfile + '.cov', 'w')
for n, record in enumerate(screed.open(dbfile)):
    if n % 10000 == 0:
        print '...', n

    print >>outfp, ">%s[cov=%d]\n%s" % (record.name, rpkms.get(record.name, 0), record.sequence)
