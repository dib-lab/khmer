#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# assert 0, "don't use this script -- use graph-size.py instead"

import khmer
import sys
import screed
import os
import subprocess
import zlib
import gzip

K = 32
HASHTABLE_SIZE = int(4e9)
THRESHOLD = 500
N_HT = 4

print 'ht size'
ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)

read_count = 0

for filename in sys.argv[1:]:

    print 'processing file: ' + filename + ' reads processed: ' + \
        str(read_count)

    for n, record in enumerate(screed.fasta.fasta_iter(open(filename))):
        seq = record['sequence']
        if len(seq) >= K:
            ht.consume(seq)

        if n % 10000 == 0:
            print '... loading', n

for filename in sys.argv[1:]:
    outfp = open(filename + '.graphsize', 'w')

    n_kept = 0

    for n, record in enumerate(screed.fasta.fasta_iter(open(filename))):
        kmer = record['sequence'][:K]
        size = ht.calc_connected_graph_size(kmer, THRESHOLD)
        if size >= THRESHOLD:
            name = record['name'].rsplit('/', 2)[0] + '.%d' % n
            print >>outfp, ">%s\n%s" % (name, record['sequence'])
            n_kept += 1

        if n % 10000 == 0:
            print '...', n, n_kept, n - n_kept + 1

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
