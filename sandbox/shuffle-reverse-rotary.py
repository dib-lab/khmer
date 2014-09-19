#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import sys
import screed
import os.path

ROTARY_SIZE = 100

prefix = os.path.basename(sys.argv[1])

fp_d = {}
for n in range(0, ROTARY_SIZE):
    num = ROTARY_SIZE - n
    fp_d[n] = open(prefix + '.%03d' % num, 'w')

total = 0
for filename in sys.argv[1:]:
    for record in screed.open(filename):
        total += 1
        if total % 10000 == 0:
            print '...', total
        loc = total % ROTARY_SIZE
        fp_d[loc].write('>%s\n%s\n' % (record.name, record.sequence))

print 'reverse-rotary shuffled %d sequences into %d files (%s.NNN)' % \
    (total, ROTARY_SIZE, prefix)
