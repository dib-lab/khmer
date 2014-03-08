#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
import sys
sys.path.insert(0, '/u/t/dev/screed')
import screed

for n, record in enumerate(screed.open(sys.argv[1])):
    if n % 10000 == 0:
        print>>sys.stderr, '...', n

    sequence = record['sequence']
    name = record['name']

    if 'N' in sequence:
        continue

    print ">" + name
    print sequence
