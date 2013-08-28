#! /usr/bin/env python
#
# This script is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import sys
from screed.fasta import fasta_iter

outfp = open(sys.argv[2], 'w')

for n, record in enumerate(fasta_iter(open(sys.argv[1]))):
    if n % 100000 == 0:
        print >>sys.stderr, '...', n

    if 'N' in record['sequence']:
        continue

    print >>outfp, '>%s\n%s' % (record['name'], record['sequence'])
