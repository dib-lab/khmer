#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import sys
import screed.fasta
import os
import khmer

K = 32
HASHTABLE_SIZE = int(8e9)
N_HT = 4

###

MAX_DENSITY = 2000

infile = sys.argv[1]
outfile = sys.argv[2]

print 'saving to:', outfile

print 'making hashtable'
ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)

print 'eating', infile
# ht.consume_fasta(infile)
ht.load(infile + '.ht')

RADIUS = 10
hist = [0.0] * 200
histcount = [0] * 200

for n, record in enumerate(screed.fasta.fasta_iter(open(infile))):
    if n % 1000 == 0:
        print '... saving', n
    seq = record['sequence']

    for pos in range(0, len(seq) - K + 1):
        density = ht.count_kmers_within_radius(seq[pos:pos + K], RADIUS,
                                               MAX_DENSITY)
        density /= float(RADIUS)
        hist[pos] += density
        histcount[pos] += 1

    if n % 1000 == 0:
        outfp = open(outfile, 'w')
        for i in range(len(hist)):
            if histcount[i]:
                print >>outfp, i, hist[i], histcount[i], hist[i] / histcount[i]
