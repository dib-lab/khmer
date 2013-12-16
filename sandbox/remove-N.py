#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import sys
import screed

outfp = open(sys.argv[2], 'w')

for n, record in enumerate(screed.open(sys.argv[1])):
    if n % 1000 == 0:
        print >>sys.stderr, '...', n

    if 'N' in record.sequence:
        continue

    print >>outfp, '>%s\n%s' % (record.name, record.sequence)
