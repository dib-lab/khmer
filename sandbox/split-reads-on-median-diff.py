#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
import sys
import math
from screed.fasta import fasta_iter
import khmer

K = 32
HASHTABLE_SIZE = int(1e9)
N_HT = 4

infile = sys.argv[1]
outfile = sys.argv[2]
outfp = open(outfile, 'w')
outfp2 = open(outfile + '.mhigh', 'w')

print 'making hashtable'
ht = khmer.new_counting_hash(K, HASHTABLE_SIZE, N_HT)

print 'eating', infile
ht.consume_fasta(infile)

print 'counting'
for n, record in enumerate(fasta_iter(open(infile))):
    if n % 10000 == 0:
        print>>sys.stderr, '...', n

    seq = record['sequence']
    if len(seq) < K:
        continue

    median, average, stddev = ht.get_median_count(seq)
    max_count = ht.get_max_count(seq)

    fp = outfp
    if max_count == 255 or median >= 128 or abs(median - average) > 60:
        fp = outfp2

    print >>fp, ">%s\n%s" % (record['name'], record['sequence'])
