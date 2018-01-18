#!/usr/bin/env python
# The following program uses khmer to
# find unique kmers to only one sequence.
import khmer
d1 = "ATGTACGGGCATTACGATTACCGATGTAG"
d2 = "ATGACCAAACTCATTACGATTAGATATAG"
ksize = 5
target_table_size = 5e5
num_tables = 4
bf = khmer.Nodetable(ksize, target_table_size, num_tables)
bf.consume(d1)
cms = khmer.Counttable(ksize, target_table_size, num_tables)
for kmer in cms.get_kmers(d2):
    if bf.get(kmer) == 0:
        cms.consume(kmer)

# If kmer is in both sequences it should not be in cms but in bf
assert cms.get('CATTA') == 0
assert bf.get('CATTA') > 0
# If kmer is in d1 but not d2 it should not be in cms but be in bf
assert cms.get('ATGTA') == 0
assert bf.get('ATGTA') > 0
# If kmer is in d2 but not d1 it should be in cms and not in bf
assert cms.get('TATAG') > 0
assert bf.get('TATAG') == 0
