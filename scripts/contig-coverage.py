import sys
import khmer

K=32
HASHTABLE_SIZE=9000
N_HT = 2

contig_file = sys.argv[1]

ht = khmer.new_counting_hash(K, HASHTABLE_SIZE, N_HT)

'''consume contigs'''
print 'consuming reads from ', contig_file

total_reads, n_consumed = ht.consume_fasta(contig_file)

print '...consumed %d reads, %d k-mers' % (total_reads, n_consumed)

ht.count('CGTGCGTTCATACTCGCAGTTTAAGAATTCTT')

print ht.get('CGTGCGTTCATACTCGCAGTTTAAGAATTCTT')
print ht.n_entries()
'''
for i in range(0, ht.n_entries()):
    n = ht.get(i)
    if n:
        print n
'''
x = ht.fasta_count_kmers_by_position(fil
 
