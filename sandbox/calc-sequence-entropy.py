#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import sys
import math
from screed.fasta import fasta_iter

for n, record in enumerate(fasta_iter(open(sys.argv[1]))):
    if n % 10000 == 0:
        print>>sys.stderr, '...', n
    seq = record['sequence']
    nA = seq.count('A')
    nC = seq.count('C')
    nG = seq.count('G')
    nT = seq.count('T')
    total = float(len(seq))

    nA /= total
    nC /= total
    nG /= total
    nT /= total

    entropy = 0
    if nA:
        entropy += - nA * math.log(nA, 2)
    if nC:
        entropy += - nC * math.log(nC, 2)
    if nG:
        entropy += - nG * math.log(nG, 2)
    if nT:
        entropy += - nT * math.log(nT, 2)

    print entropy
