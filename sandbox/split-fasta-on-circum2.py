#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
# saved for posterity; DO NOT USE.
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

infile = sys.argv[1]
outprefix = sys.argv[2]

lowfile = outprefix + '.low'
highfile = outprefix + '.high'
circfile = outprefix + '.circ'

print 'saving low-density to:', lowfile
print 'saving high-density to:', highfile
print 'saving all densities to:', circfile

print 'making hashtable'
ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)

lowfp = open(lowfile, 'w')
highfp = open(highfile, 'w')
circfp = open(circfile, 'w')

print 'eating', infile
ht.consume_fasta(infile)

start = RADIUS
incr = 2 * RADIUS

for n, record in enumerate(screed.fasta.fasta_iter(open(infile),
                                                   parse_description=False)):
    if n % 10000 == 0:
        print '... saving', n

    seq = record['sequence']

    # between [RADIUS:-RADIUS] kmers, calculate circumference every 2*RADIUS
    end = len(seq) - K + 1 - incr / 2

    is_high = False
    circums = []
    for pos in range(end, start, -incr):
        circum = ht.count_kmers_on_radius(seq[pos:pos + K], RADIUS, MAX_VOLUME)

        circums.append((circum, pos))
        print >>circfp, circum

        if circum >= MAX_CIRCUM:
            is_high = True

    # did we find a circumference >= our max?  If so, try to trim the reads.
    if is_high:
        # trim from the back:
        for (circum, pos) in circums:
            if circum < MAX_CIRCUM:
                break

        chop_end = pos

        # trim from the front:
        circums.reverse()
        for (circum, pos) in circums:
            if circum < MAX_CIRCUM:
                break
        chop_start = pos

        # make sure we're not missing anything in the middle that wasn't
        # trimmed on either side:
        if max([circ for (circ, _) in circums]) < MAX_CIRCUM and \
                chop_end - chop_start >= 0:

            # do the trimming & rewrite the name
            sequence = record['sequence'][chop_start:chop_end + K]
            record['name'] = record['name'] + '\tTRUNC:%d-%d' % (chop_start,
                                                                 chop_end + K)

            # clear the 'is_high' flag if there are no more high-circum
            # k-mers left.
            is_high = False

    # sort "high circumference" and "low" circumference sequences separately.
    if is_high:
        fp = highfp
    else:
        fp = lowfp

    print >>fp, '>%s\n%s' % (record['name'], record['sequence'])

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
