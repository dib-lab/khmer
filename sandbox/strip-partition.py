#! /usr/bin/env python
import screed
import sys

for record in screed.open(sys.argv[1]):
    name = record['name']
    sequence = record['sequence']

    name = name.split()[0]
    
    print '>%s\n%s' % (name, sequence,)
