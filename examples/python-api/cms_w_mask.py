#!/usr/bin/env python

# A demonstration of using khmer to populate a count-min sketch (cms)
# with a mask.
# Put all the kmers from dataset2 into cms except those that
# are shared with dataset1.
# Typically khmer accrues a small false positive rate in order to save
# substantially on memory requirements.

import khmer

ksize = 6
target_table_size = 5e8
num_tables = 4

# strings of your datasets
dataset1 = 'GCTGCACCGATGTACGCAAAGCTATTTAAAACCATAACTATTCTCACTTA'
dataset2 = 'CCTGCACCGACGTACGCTATGCTATTGAAGACCATTAGTAGGCTCACTCC'

# create a bloom filter
bloomfilter = khmer.Nodetable(ksize, target_table_size, num_tables)
# load dataset1 and store all the kmers
bloomfilter.consume(dataset1)

cms = khmer.Counttable(ksize, target_table_size, num_tables)

# for every kmer in dataset2
for kmer in cms.get_kmers(dataset2):
    if bloomfilter.get(kmer) == 0:  # kmers unique to cms
        cms.consume(kmer)

# this kmer is in dataset2 (cms), but not dataset1
assert cms.get('CCTGCA') > 0

# this kmer is in dataset1 (bloomfilter), but not dataset2
assert bloomfilter.get('GCTGCA') > 0

# this kmer is in both datasets, should not be in cms
assert cms.get('GTACGC') == 0
