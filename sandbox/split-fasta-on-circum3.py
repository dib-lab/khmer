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

###

RADIUS = 2
MAX_CIRCUM = 4                            # 4 seems to eliminate lump in 1m.fa
MAX_VOLUME = 200

repfile = sys.argv[1]
infile = sys.argv[2]
outprefix = sys.argv[3]

lowfile = outprefix + '.low'
highfile = outprefix + '.high'

print 'saving low-density to:', lowfile
print 'saving high-density to:', highfile

print 'making hashtable'
ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)

lowfp = open(lowfile, 'w')
highfp = open(highfile, 'w')

print 'eating', infile
ht.consume_fasta(repfile)

incr = 2 * RADIUS

for n, record in enumerate(screed.fasta.fasta_iter(open(infile),
                                                   parse_description=False)):
    if n % 10000 == 0:
        print '... saving', n

    name = record['name']
    seq = record['sequence']

    # calculate circumference for every point.
    end = len(seq) - K
    is_high = False

    pos = 0
    for pos in range(0, end, incr):
        circum = ht.count_kmers_on_radius(seq[pos:pos + K], RADIUS, MAX_VOLUME)

        if circum >= MAX_CIRCUM:
            is_high = True
            break

    # ok. sequence has high-radius k-mers; can we trim them off?
    if is_high and pos >= incr:
        pos -= incr

        # find last k-mer with a low radius:
        i = 1
        for i in range(1, incr):
            circum = ht.count_kmers_on_radius(seq[pos + i:pos + i + K],
                                              RADIUS, MAX_VOLUME)
            if circum >= MAX_CIRCUM:
                break

        pos += i - 1

        # now trim sequence:
        seq = seq[:pos + K]
        is_high = False
        name += "\tTRUNC.%d" % pos

    # sort "high circumference" and "low" circumference sequences separately.
    if is_high:
        fp = highfp
    else:
        fp = lowfp

    print >>fp, '>%s\n%s' % (name, seq)
