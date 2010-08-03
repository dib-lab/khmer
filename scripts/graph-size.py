import khmer
import sys

K = 32
HASHTABLE_SIZE=4**12+1

infile = sys.argv[1]
threshold = int(sys.argv[2])
outfile = sys.argv[3]

print 'creating ht'
ht = khmer.new_hashtable(K, HASHTABLE_SIZE)
print 'eating fa', infile
total_reads, n_consumed = ht.consume_fasta(infile)

print 'hashtable occupancy:', ht.n_occupied() / float(HASHTABLE_SIZE)
#fp = open('aaa.1', 'w')
#total = 0
#for n, i in enumerate(ht.graphsize_distribution(10000)):
#    total += i
#    print >>fp, n, i, total
#fp.close()

print 'trimming to', threshold
ht.trim_graphs(infile, threshold, outfile)

#print 'hashtable occupancy:', ht.n_occupied() / float(HASHTABLE_SIZE)

#fp = open('aaa.2', 'w')
#total = 0
#for n, i in enumerate(ht.graphsize_distribution(5000)):
#    total += i
#    print >>fp, n, i, total
#fp.close()
