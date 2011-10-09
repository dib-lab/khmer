#!/usr/bin/env python

import sys
import screed
from screed import fasta

filein = sys.argv[1]

fp = open(filein)

lengths = [0] * 100
for n, record in enumerate(fasta.fasta_iter(fp)):
    length = len(record['sequence']) - 32
    lengths[length] += 1

for n, i in enumerate(lengths):
    print n, i
