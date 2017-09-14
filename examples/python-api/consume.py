#!/usr/bin/env python

# A demonstration of khmer's primary sequence loading function.

import khmer
import sys

ksize = 21
target_table_size = 5e8
num_tables = 4

counts = khmer.Counttable(ksize, target_table_size, num_tables)
nseqs, nkmers = counts.consume_seqfile(sys.argv[1])
print('Loaded', nseqs, 'sequences and', nkmers, 'k-mers from', sys.argv[1])

print('The kmer "CAGCGCCGTGTTGTTGCAATT" appears',
      counts.get('CAGCGCCGTGTTGTTGCAATT'), 'times in the data')
print('The kmer "GATTACAGATTACAGATTACA" appears',
      counts.get('GATTACAGATTACAGATTACA'), 'times in the data')
