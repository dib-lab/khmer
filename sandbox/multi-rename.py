#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import screed
import sys

CUTOFF = 200


n = 0
for filename in sys.argv[1:]:
    for record in screed.open(filename):
        if len(record.sequence) >= CUTOFF:
            n += 1
            print '>%s %s\n%s' % (n, record.name, record.sequence)
