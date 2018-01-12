########################################
#The following program uses khmer to 
#find unique kmers to only one sequence. 
########################################
#!/usr/bin/env python
import khmer
d1 = "ATGTACGGGCATTACGATTACCGATGTAG"
d2 = "ATGACCAAACTCATTACGATTAGATATAG"
ksize = 5
target_table_size = 5e8
num_tables = 4
bf = khmer.Nodetable(ksize, target_table_size, num_tables)
bf.consume(d1)
cms = khmer.Counttable(ksize, target_table_size, num_tables)
for kmer in cms.get_kmers(d2):
    if bf.get(kmer) == 0:
        cms.consume(kmer)
assert cms.get('CATTA') == 0
