import khmer
import sys

K = 17
HASHTABLE_SIZE=4**K

infile = sys.argv[1]
threshold = int(sys.argv[2])
outfile = sys.argv[3]

print 'creating ht'
ht = khmer.new_hashtable(K, HASHTABLE_SIZE)
print 'eating fa', infile
total_reads, n_consumed = ht.consume_fasta(infile)

print 'hashtable occupancy:', ht.n_occupied() / float(HASHTABLE_SIZE)
fp = open('aaa.1', 'w')
for n, i in enumerate(ht.graphsize_distribution(500)):
    print >>fp, n, i
ht.clear_marks()
    
print 'trimming to', threshold
ht.trim_graphs(threshold)
ht.clear_marks()
print 'hashtable occupancy:', ht.n_occupied() / float(HASHTABLE_SIZE)

fp = open('aaa.2', 'w')
for n, i in enumerate(ht.graphsize_distribution(500)):
    print >>fp, n, i
    
print 'filtering'
minmax = ht.fasta_file_to_minmax(infile, total_reads)
readmask = ht.filter_fasta_file_any(minmax, 1)

readmask.filter_fasta_file(infile, outfile)
