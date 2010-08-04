import khmer
import sys
import screed

K = 32
HASHTABLE_SIZE=4**12+1

infile = sys.argv[1]
threshold = int(sys.argv[2])
outfile = sys.argv[3]

print 'creating ht'
ht = khmer.new_hashtable(K, HASHTABLE_SIZE)
print 'eating fa', infile
total_reads, n_consumed = ht.consume_fasta(infile)
outfp = open(outfile, 'w')

#print 'hashtable occupancy:', ht.n_occupied() / float(HASHTABLE_SIZE)

for record in screed.fasta.fasta_iter(open(infile)):
    kmer = record['sequence'][:K]
    size = ht.calc_connected_graph_size(kmer, threshold)
    if size >= threshold:
        print >>outfp, ">%s\n%s" % (record['name'], record['sequence'])

#print 'trimming to', threshold
#ht.trim_graphs(infile, threshold, outfile)

#print 'hashtable occupancy:', ht.n_occupied() / float(HASHTABLE_SIZE)

#fp = open('aaa.2', 'w')
#total = 0
#for n, i in enumerate(ht.graphsize_distribution(5000)):
#    total += i
#    print >>fp, n, i, total
#fp.close()
