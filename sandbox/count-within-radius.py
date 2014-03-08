#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import sys
import screed.fasta
import os
import khmer

K = 32
HASHTABLE_SIZE = int(8e9)
N_HT = 4
RADIUS = 100

###

MAX_DENSITY = 2000

infile = sys.argv[1]
outfile = sys.argv[2]
if len(sys.argv) > 3:
    RADIUS = int(sys.argv[3])

print 'saving to:', outfile

print 'making hashtable'
ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)

print 'eating', infile
ht.consume_fasta(infile)

print 'loading'
ht.save(outfile + '.ht')

outfp = open(outfile, 'w')
for n, record in enumerate(screed.open(infile)):
    if n % 10000 == 0:
        print '... saving', n
    seq = record['sequence']

    for pos in range(0, len(seq), 200):
        subseq = seq[pos:pos + 200]

        middle = (len(subseq) - K + 1) / 2

        density = ht.count_kmers_within_radius(
            subseq[middle:middle + K], RADIUS,
            MAX_DENSITY)
        density /= float(RADIUS)

        print >>outfp, '>%s d=%.3f\n%s' % (record['name'], density, subseq)
