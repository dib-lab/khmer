#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
## using a HyperLogLog counter and a bloom filter to count unique kmers,
## comparing results

import khmer
import sys
from screed.fasta import fasta_iter


filename = sys.argv[1]
K = int(sys.argv[2])  # size of kmer

ERROR_RATE = .01

hllcpp = khmer.new_hll_counter(ERROR_RATE, K)

hllcpp.consume_fasta(filename)

cpp_estimate = hllcpp.estimate_cardinality()

print(cpp_estimate)
