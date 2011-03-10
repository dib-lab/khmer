#! /usr/bin/env python
import sys
import screed, khmer

K = 32

infile = sys.argv[1]

ht = khmer.new_hashbits(K, 1, 1)
ht.consume_partitioned_fasta(infile)

for n, record in enumerate(screed.open(infile)):
    if n % 10000 == 0:
        print '... checking', n
    assert ht.is_single_partition(record.sequence)
