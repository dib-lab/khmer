#! /usr/bin/env python
#
# This script is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import khmer
import sys
import screed

K = 32
HASHTABLE_SIZE = int(1e9)


def get_neighbors(kmer):
    neighbors = []
    bases = ['A', 'C', 'G', 'T']

    for n, base in enumerate(kmer):
        for new_base in bases:
            if new_base == base:
                continue

            new_kmer = kmer[0:n] + new_base + kmer[n + 1:]
            neighbors.append(new_kmer)

    return neighbors

ht = khmer.new_hashtable(K, HASHTABLE_SIZE)

filename = sys.argv[1]

ht.consume_fasta(filename)

for n, record in enumerate(screed.fasta.fasta_iter(open(filename))):
    read = record['sequence']
    name = record['name']

    if len(read) < K:
        continue

    seq_len = len(read)
    for n in range(0, seq_len + 1 - K):
        kmer = read[n:n + K]
        kmer_count = ht.get(kmer)

        if kmer_count == 1:
            neighbors = get_neighbors(kmer)

            for neighbor in neighbors:
                neig_count = ht.get(neighbor)
                if neig_count > 1:
                    read = read[0:n] + neighbor + read[n + K:]

                    break

    print ">" + name
    print read
