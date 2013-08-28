#! /usr/bin/env python
#
# This script is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import sys
import screed
import khmer

K = 32

infile = sys.argv[1]

ht = khmer.new_hashbits(K, 1, 1)
ht.consume_partitioned_fasta(infile)

for n, record in enumerate(screed.open(infile)):
    if n % 10000 == 0:
        print '... checking', n
    assert ht.is_single_partition(record.sequence)
