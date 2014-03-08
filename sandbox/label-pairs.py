#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
import screed
import sys
import itertools

s1_file = sys.argv[1]

n = 0
for r in screed.open(s1_file):
    n += 1

    name = r.name
    if n % 2 == 0:
        if not name.endswith('/2'):
            name += '/2'
    elif not name.endswith('/1'):
        name += '/1'

    print '>%s\n%s' % (name, r.sequence)
