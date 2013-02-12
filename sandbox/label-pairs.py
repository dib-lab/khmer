#! /usr/bin/env python
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
