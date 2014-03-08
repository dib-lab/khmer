#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
import khmer
import sys
import screed

K = 20
HASHTABLE_SIZE = int(2.5e8)
N_HT = 4

UNIQUE_LEN = 100
UNIQUE_F = 0.9

OUTPUT_WINDOW = 100
OUTPUT_OVERLAP = 10

kh = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)

discarded = 0
kept_kmers = 0
total_kmers = 0

total_out = 0
for filename in sys.argv[1:]:
    n_out = 0
    for n, record in enumerate(screed.open(filename)):
        if n > 0 and n % 10000 == 0:
            print >>sys.stderr, '...', n, discarded
            print >>sys.stderr, '==>', total_kmers, kept_kmers, int(
                float(kept_kmers) / float(total_kmers) * 100.)
        seq = record.sequence
        seq = seq.replace('N', 'G')

        paths = kh.extract_unique_paths(seq, UNIQUE_LEN, UNIQUE_F)

        kh.consume(seq)
        total_kmers += len(seq) - K + 1

        if not len(paths):
            discarded += 1
            continue

        for i, path in enumerate(paths):
            n_out += 1

            if len(path) < OUTPUT_WINDOW:
                total_out += 1
                print '>%d\n%s' % (total_out, path)
                continue

            for start in range(0, len(path) - OUTPUT_WINDOW + 1,
                               OUTPUT_OVERLAP):
                total_out += 1
                subpath = path[start:start + OUTPUT_WINDOW]
                print '>%d\n%s' % (total_out, subpath)

    print >>sys.stderr, '%d for %s' % (n_out, filename)
