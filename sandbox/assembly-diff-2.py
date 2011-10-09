#! /usr/bin/env python
import sys, khmer, screed, os

K=20
HASHTABLE_SIZE=int(4e9)
N_HT=4

UNIQUE_LEN=100
UNIQUE_F=0.9

filename1 = sys.argv[1]
filename2 = sys.argv[2]
uniq2 = open(os.path.basename(sys.argv[2]) + '.uniq', 'w')

kh = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)
for n, record in enumerate(screed.open(filename1)):
    if n % 10000 == 0:
        print '...', filename1, n
    seq = record.sequence.upper().replace('N', 'G')
    kh.consume(seq)

path_n = 0
for n, record in enumerate(screed.open(filename2)):
    if n % 10000 == 0:
        print '...', filename2, n
    seq = record.sequence.upper().replace('N', 'G')
    paths = kh.extract_unique_paths(seq, UNIQUE_LEN, UNIQUE_F)
    kh.consume(seq)

    for path in paths:
        path_n += 1
        print >>uniq2, '>%s\n%s' % (path_n, path)
