#! /usr/bin/env python
import sys
import screed

for record in screed.open(sys.argv[1]):
    print '>%s\n%s' % (record.name, record.sequence)
