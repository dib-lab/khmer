#! /usr/bin/env python
#
# This script is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import sys
sys.path.insert(0, '/u/t/dev/screed')
import screed
from screed.fastq import fastq_iter

import khmer

hasher = khmer.HashtableIntersect(15, *khmer.PRIMES_1m)

filename = '/scratch/gpgc/iowa/850.2.fq'

for n, record in enumerate(fastq_iter(open(filename))):
    if n % 10000 == 0:
        print>>sys.stderr, '...', n
    if n > 10 ** 6:
        break

    sequence = record['sequence']
    if 'N' in sequence:
        continue

    print '>%s\n%s' % (record['name'], record['sequence'])
