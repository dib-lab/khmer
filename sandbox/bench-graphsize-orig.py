#! /usr/bin/env python
#
# This script is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import khmer
import sys
import screed
import os
import subprocess
import zlib
import gzip

K = 14
HASHTABLE_SIZE = int(1e9)
THRESHOLD = 100

print 'ht size'
ht = khmer.new_hashbits(K, HASHTABLE_SIZE, 1)

read_count = 0

for filename in sys.argv[1:]:
    ht.consume_fasta(filename)

    if 0:
        print 'processing file: ' + filename + ' reads processed: ' + \
            str(read_count)
        for n, record in enumerate(screed.fasta.fasta_iter(open(filename))):
            read_count += 1
            ht.consume(record['sequence'])

            if n % 10000 == 0:
                print '...', n

for filename in sys.argv[1:]:
    outfp = open(filename[:-3] + '.graphsize', 'w')

    n_kept = 0

    for n, record in enumerate(screed.fasta.fasta_iter(open(filename))):
        kmer = record['sequence'][:K]
        size = ht.calc_connected_graph_size(kmer, THRESHOLD)
        if size >= THRESHOLD:
            print >>outfp, ">%s\n%s" % (record['name'], record['sequence'])
            n_kept += 1

        if n % 10000 == 0:
            print '...', n, n_kept, n - n_kept + 1

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
