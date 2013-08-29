#! /usr/bin/env python
import sys
import khmer
import screed
import os

K = 20
HASHTABLE_SIZE = int(2.5e8)
N_HT = 4

THRESHOLD = 0.9

kh = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)
kh.consume_fasta(sys.argv[1])

outfp = open(os.path.basename(sys.argv[1]) + '.vector', 'w')

for record in screed.open(sys.argv[2]):
    n = 0
    n_present = 0

    path = record.sequence

    n = len(path) - K + 1
    for i in range(n):
        if kh.get(path[i:i + K]):
            n_present += 1

    if n_present / float(n) >= THRESHOLD:
        print >>outfp, '1',
    else:
        print >>outfp, '0',
