#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import sys
import khmer

files = sys.argv[2:]

total_reads = len(files) * [0]
n_consumed = len(files) * [0]
n_seq_kept = len(files) * [0]

print 'loading ht'
ht = khmer.new_counting_hash(1, 1, 1)

ht.load(sys.argv[1])

for i, infile in enumerate(files):
   print 'outputting', infile + '.freq'
   ht.output_fasta_kmer_pos_freq(infile, infile + ".freq")
