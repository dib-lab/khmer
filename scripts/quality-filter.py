#!/usr/bin/env python

import sys
import screed
from screed import fastq

# python quality-filter.py <input fastq file> <output filtered fastq file>
# MINLENGTH is the minimum lenth of read desired.  NCALLS is the percentage of a read with 'N' base calls for which if read has greater, it will be removed. 

MINLENGTH = 30
NCALLS = 0.70


filein = sys.argv[1]
fileout = sys.argv[2]

fp = open(filein)
fw = open(fileout, 'w')

count=0
for n, record in enumerate(fastq.fastq_iter(fp)):
    name = record['name']
    sequence = record['sequence']
    accuracy = record['accuracy']
    
    if float(sequence.count('N')) / len(sequence) > NCALLS:
        continue

    else:
        trim = accuracy.find('B')

        if trim > MINLENGTH or (trim == -1 and len(sequence) > MINLENGTH):
            if trim == -1:
                fw.write('@%s\n%s\n+\n%s\n' % (name, sequence, accuracy))
            else:
                fw.write('@%s\n%s\n+\n%s\n' % (name, sequence[:trim], accuracy[:trim]))
            count += 1

    if n % 1000 == 0:
        print 'scanning', n

print 'Original Number of Reads', n + 1
print 'Final Number of Reads', count
print 'Total Filtered', n + 1  - int(count)
