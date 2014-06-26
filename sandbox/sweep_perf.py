#! /usr/bin/env python2

import khmer
import screed
import sys
import time

R = int(sys.argv[1])
print R
K = 20
test_file = '../tests/test-data/biglump-random-20-a.fa'

ht = khmer.new_hashbits(K, 1e9, 4)
ht.consume_fasta_and_tag_with_colors(test_file)

N = 10
for n, record in enumerate(screed.open(test_file)):
    if n > N:
        break
    print '*' * 40
    print '{} k-mers in sequence'.format(len(record.sequence) - K + 1)

    stime = time.clock()
    colors = ht.sweep_color_neighborhood(record.sequence, R)
    etime = time.clock()

    print 'traversal took {} seconds'.format(etime - stime)
    print 'found {} colors'.format(len(colors))
