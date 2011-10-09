#! /usr/bin/env python
import sys
import math
from screed.fasta import fasta_iter
import khmer

K = 32
HASHTABLE_SIZE=int(1e9)
N_HT = 4

infile = sys.argv[1]
outfile = sys.argv[2]
outfp = open(outfile, 'w')

print 'making hashtable'
ht = khmer.new_counting_hash(K, HASHTABLE_SIZE, N_HT)

print 'eating', infile
ht.consume_fasta(infile)

print 'counting'
for n, record in enumerate(fasta_iter(open(infile))):
    if n % 10000 == 0:
        print>>sys.stderr, '...', n
        
    seq = record['sequence']
    if len(seq) < K:
        continue

    x = []
    for pos in range(0, len(seq) - K + 1):
        x.append(ht.get(seq[pos:pos+K]))

    print >>outfp, '>%s\n%s' % (record['name'], record['sequence'])
    print >>outfp, " ".join(map(str, x))
    
    median, average, stddev = ht.get_median_count(seq)
    max_count = ht.get_max_count(seq)
    min_count = ht.get_min_count(seq)

    print >>outfp, "%d %.2f %.2f %d %d" % (median, average, stddev,
                                           min_count, max_count)
