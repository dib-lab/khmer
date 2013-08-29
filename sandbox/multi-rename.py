#! /usr/bin/env python
import screed
import sys

CUTOFF = 200


n = 0
for filename in sys.argv[1:]:
    for record in screed.open(filename):
        if len(record.sequence) >= CUTOFF:
            n += 1
            print '>%s %s\n%s' % (n, record.name, record.sequence)
