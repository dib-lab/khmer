#! /usr/bin/env python
import screed
import sys
import os.path

fp1 = open(os.path.basename(sys.argv[1]) + '.1', 'w')
fp2 = open(os.path.basename(sys.argv[1]) + '.2', 'w')

for n, record in enumerate(screed.open(sys.argv[1])):
    if n % 10000 == 0:
        print >>sys.stderr, '...', n

    name = record.name
    if name.endswith('/1'):
        print >>fp1, '>%s\n%s' % (record.name, record.sequence,)
    elif name.endswith('/2'):
        print >>fp2, '>%s\n%s' % (record.name, record.sequence,)
