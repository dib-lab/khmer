#! /usr/bin/env python
import screed, sys

CUTOFF=200


n = 0
for filename in sys.argv[1:]:
   for record in screed.open(filename):
       n += 1
       print '>%s %s\n%s' % (n, record.name, record.sequence)
