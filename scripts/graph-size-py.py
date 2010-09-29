import khmer
import sys
import screed

K = 32
HASHTABLE_SIZE=int(256e9)
THRESHOLD=100

print 'creating ht'
ht = khmer.new_hashbits(K, HASHTABLE_SIZE, 2)

for filename in sys.argv[1:]:
    print 'eating fa', filename
    total_reads, n_consumed = ht.consume_fasta(filename)

#ht.save(sys.argv[1] + '.graphsize.ht')

for filename in sys.argv[1:]:
    outfp = open(filename + '.graphsize', 'w')

    n_kept = 0
    for n, record in enumerate(screed.fastq.fastq_iter(open(filename))):
        kmer = record['sequence'][:K]
        size = ht.calc_connected_graph_size(kmer, THRESHOLD)
        if size >= THRESHOLD:
            print >>outfp, ">%s\n%s" % (record['name'], record['sequence'])
            n_kept += 1

        if n % 10000 == 0:
            print '...', n, n_kept, n - n_kept

    
