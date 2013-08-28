#! /usr/bin/env python
#
# This script is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import sys
import screed.fasta
import os
import khmer

K = 31                                  # use K-1 for assembly K
HASHTABLE_SIZE = int(4e9)
N_HT = 4

###

MAX_DEGREE = int(sys.argv[4])

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

print 'eating', repfile
ht.consume_fasta(repfile)

for n, record in enumerate(screed.fasta.fasta_iter(open(infile),
                                                   parse_description=False)):
    if n % 10000 == 0:
        print '... saving', n

    name = record['name']
    seq = record['sequence']
    trim_seq, trim_at = ht.trim_on_degree(seq, MAX_DEGREE)

    # sort high degree & low degree sequences separately; ignore trimmed
    # component.

    if trim_at > K:
        print >>lowfp, '>%s %d %d %d\n%s' % (
            name, trim_at, len(seq), len(trim_seq), trim_seq)
    else:
        print >>highfp, '>%s\n%s' % (name, seq)
