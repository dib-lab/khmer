#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
import sys
import khmer

filename = sys.argv[1]
k = int(sys.argv[2])

ht = khmer.new_hashtable(k, 4 ** k)
ht.consume_fasta(filename)

N = ht.n_occupied()

for i in range(0, 4):
    print 'q', i, ht.n_occupied(i * 4 ** (k - 1), (i + 1) * 4 ** (k - 1))
