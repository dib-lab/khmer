#! /usr/bin/env python
import sys, screed

K = 32

for filename in sys.argv[1:]:
   n_kmers = 0
   n_reads = 0

   for record in screed.open(filename):
      n_kmers += len(record.sequence) - K + 1
      n_reads += 1

   fp = open(filename + '.stats', 'w')
   fp.write('K=%d\nn_kmers: %d\nn_reads: %d\n' % (K, n_kmers, n_reads))
   fp.close()
