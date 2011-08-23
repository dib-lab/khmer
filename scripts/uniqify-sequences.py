#! /usr/bin/env python
import khmer, sys, screed

K=20
HASHTABLE_SIZE=int(1e9)
N_HT=4

UNIQUE_LEN=100
UNIQUE_F=0.9

kh = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)

discarded = 0
kept_kmers = 0
total_kmers = 0
for n, record in enumerate(screed.open(sys.argv[1])):
    if n > 0 and n % 10000 == 0:
        print >>sys.stderr, '...', n, discarded
        print >>sys.stderr, '==>', total_kmers, kept_kmers, int(float(kept_kmers) / float(total_kmers) * 100.)
    seq = record.sequence
    seq = seq.replace('N', 'G')

    paths = kh.extract_unique_paths(seq, UNIQUE_LEN, UNIQUE_F)

    kh.consume(seq)
    total_kmers += len(seq) - K + 1

    if not len(paths):
        discarded += 1
        continue

    for i, path in enumerate(paths):
        kept_kmers += len(path) - K + 1
        print '>%s.%i\n%s' % (record.name, i, path)
