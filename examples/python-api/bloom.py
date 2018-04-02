#!/usr/bin/env python

# A demonstration of using khmer to query a dataset for a k-mer. Typically
# khmer accrues a small false positive rate in order to save substantially on
# memory requirements.

import khmer

ksize = 21
target_table_size = 5e8
num_tables = 4

bloomfilter = khmer.Nodetable(ksize, target_table_size, num_tables)
bloomfilter.consume('GCTGCACCGATGTACGCAAAGCTATTTAAAACCATAACTATTCTCACTTA')

print('count for "GCTGCACCGATGTACGCAAAG" is',
      bloomfilter.get('GCTGCACCGATGTACGCAAAG'))

bloomfilter.count('GCTGCACCGATGTACGCAAAG')

print('count for "GCTGCACCGATGTACGCAAAG" is',
      bloomfilter.get('GCTGCACCGATGTACGCAAAG'))

print('count for "GATTACAGATTACAGATTACA" is',
      bloomfilter.get('GATTACAGATTACAGATTACA'))
