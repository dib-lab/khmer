#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import sys
import khmer
import screed

K = 32

ht = khmer.load_counting_hash(sys.argv[1])

for record in screed.open(sys.argv[2]):
    for pos in range(len(record.sequence) - K + 1):
        if ht.get(record.sequence[pos:pos + K]) > 40000:
            print record.sequence[pos:pos + K]
            sys.exit(0)
