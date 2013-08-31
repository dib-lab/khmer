#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#

import sys
import screed
import gzip

# python quality-trim.py <input fastq file> <output filtered fastq file>
# MINLENGTH is the minimum length of read desired.  NCALLS is the
# percentage of a read with 'N' base calls for which if read has greater,
# it will be removed.

MINLENGTH = 30

filein = sys.argv[1]
fileout = sys.argv[2]

fw = gzip.open(fileout, 'w')

count = 0
for n, record in enumerate(screed.open(filein)):
    name = record['name']
    sequence = record['sequence']
    accuracy = record['accuracy']

    sequence = sequence.rstrip('N')
    accuracy = accuracy[:len(sequence)]

    if 'N' in sequence:
        continue
    else:
        trim = accuracy.find('B')

        if trim > MINLENGTH or (trim == -1 and len(sequence) > MINLENGTH):
            if trim == -1:
                fw.write('@%s\n%s\n+\n%s\n' % (name, sequence, accuracy))
            else:
                fw.write('@%s\n%s\n+\n%s\n' % (name,
                         sequence[:trim], accuracy[:trim]))
            count += 1

    if n % 1000 == 0:
        print 'scanning', n

print 'Original Number of Reads', n + 1
print 'Final Number of Reads', count
print 'Total Filtered', n + 1 - int(count)
