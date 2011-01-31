#! /usr/bin/env python
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
