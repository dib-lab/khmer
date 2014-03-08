#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
import khmer
import screed
import sys

K = 32

ht = khmer.new_hashbits(K, 1000000000, 8)

unique_count = 0

filename = sys.argv[1]

for n, record in enumerate(screed.fasta.fasta_iter(open(filename))):
   read = record['sequence']
   name = record['name']

   if len(read) < K:
      continue

   seq_len = len(read)
   for n in range(0,seq_len+1-K):
      kmer = read[n:n+K]

      if ht.get(kmer) == 0:
         unique_count += 1
         ht.count(kmer)

print unique_count
print ht.n_occupied()
