#! /usr/bin/env python
#
# This script is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
# this only works for small files! it loads everything into mem.
import sys
from screed.fasta import fasta_iter
import random

d = dict([(r['name'], r['sequence']) for r in fasta_iter(open(sys.argv[1]))])

ks = d.keys()
random.shuffle(ks)

for k in ks:
    s = d[k]

    print '>%s\n%s' % (k, s)
