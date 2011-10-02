#! /usr/bin/env python
import sys, khmer, screed, os

K=20
HASHTABLE_SIZE=int(2.5e8)
N_HT=4

UNIQUE_LEN=100
UNIQUE_F=0.9

filename1 = sys.argv[1]
filename2 = sys.argv[2]
uniq2 = open(os.path.basename(sys.argv[2]) + '.uniq', 'w')

kh = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)
kh.consume_fasta(filename1)

n = 0
for record in screed.open(filename2):
    seq = record.sequence
    paths = kh.extract_unique_paths(seq, UNIQUE_LEN, UNIQUE_F)
    kh.consume(seq)

    for path in paths:
        n += 1
        print >>uniq2, '>%s\n%s' % (n, path)
