#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import sys
import khmer
from screed.fasta import fasta_iter

K = 32
HASHTABLE_SIZE = int(8e8)
N_HASHTABLES = 4

total_kmers = 0

ht = khmer.new_hashbits(32, HASHTABLE_SIZE, N_HASHTABLES)
pidset = set()
for n, record in enumerate(
        fasta_iter(open(sys.argv[1]), parse_description=False)):
    pid = record['name'].rsplit('\t', 1)[1]
    pidset.add(pid)

    ht.consume(record['sequence'])
    total_kmers += len(record['sequence']) - K + 1
unique_kmers = ht.n_unique_kmers()

print 'n partitions:', len(pidset)
print 'unique kmers:', unique_kmers
print 'total kmers:', total_kmers
print 'average kmer coverage: %.2f' % (total_kmers / float(unique_kmers))

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
