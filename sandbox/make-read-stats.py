#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
import sys
import screed

K = 32

for filename in sys.argv[1:]:
    n_kmers = 0
    n_reads = 0

    for record in screed.open(filename):
        n_kmers += len(record.sequence) - K + 1
        n_reads += 1

    fp = open(filename + '.stats', 'w')
    fp.write('K=%d\nn_kmers: %d\nn_reads: %d\n' % (K, n_kmers, n_reads))
    fp.close()
