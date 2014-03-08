#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
import sys
import khmer
from screed.fasta import fasta_iter

K = 32
HASHTABLE_SIZE = int(1e9)
N_HASHTABLES = 4

ht = khmer.new_counting_hash(32, HASHTABLE_SIZE, N_HASHTABLES)

filename = sys.argv[1]

ht.consume_fasta(filename)
total, count, mean = ht.get_kmer_abund_mean(filename)
abs_dev = ht.get_kmer_abund_abs_deviation(filename, mean)

print total, count, mean, abs_dev

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
