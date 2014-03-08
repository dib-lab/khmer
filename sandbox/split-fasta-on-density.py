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

RADIUS = 10
MAX_DENSITY = 31

infile = sys.argv[1]
outprefix = sys.argv[2]

lowfile = outprefix + '.low'
highfile = outprefix + '.high'
densfile = outprefix + '.dens'

print 'saving low-density to:', lowfile
print 'saving high-density to:', highfile
print 'saving all densities to:', densfile

print 'making hashtable'
ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)

lowfp = open(lowfile, 'w')
highfp = open(highfile, 'w')
densfp = open(densfile, 'w')

print 'eating', infile
ht.consume_fasta(infile)

start = RADIUS / 2
incr = RADIUS

for n, record in enumerate(screed.fasta.fasta_iter(open(infile),
                                                   parse_description=False)):
    if n % 10000 == 0:
        print '... saving', n

    seq = record['sequence']
    end = len(seq) - K + 1 - incr / 2

    is_high = False
    densities = []
    for pos in range(end, start, -incr):
        density = ht.count_kmers_within_radius(seq[pos:pos + K], RADIUS,
                                               MAX_DENSITY)

        densities.append((density, pos))
        print >>densfp, density

        if density >= MAX_DENSITY:
            is_high = True

    if is_high:
        for (density, pos) in densities:
            if density < MAX_DENSITY:
                break

        chop_end = pos + K

        densities.reverse()
        for (density, pos) in densities:
            if density < MAX_DENSITY:
                break
        chop_start = pos

        if chop_end - chop_start > 0:
            sequence = record['sequence'][chop_start:chop_end]
            record['name'] = record['name'] + '\tTRUNC:%d-%d' % (chop_start,
                                                                 chop_end)
            is_high = False

    if is_high:
        fp = highfp
    else:
        fp = lowfp

    print >>fp, '>%s\n%s' % (record['name'], record['sequence'])
