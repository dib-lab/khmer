import khmer
import sys
import screed
import os
import subprocess
import zlib
import gzip

K = 30
HASHTABLE_SIZE=4**17+1
THRESHOLD=100

print 'ht size'
ht = khmer.new_hashbits(K, HASHTABLE_SIZE, 1)

read_count = 0

fd = open(sys.argv[1], 'w')

for filename in sys.argv[2:]:

    print 'processing file: ' + filename + ' reads processed: ' + str(read_count)

    for n, record in enumerate(screed.fastq.fastq_iter(gzip.open(filename))):

        read_count += 1
        ht.consume(record['sequence'])
    
        if read_count % 100000 == 0:
            fd.write(str(read_count) + " " + str(ht.n_occupied()) + " " + str(ht.n_occupied() / float(HASHTABLE_SIZE)) +'\n')
            fd.flush()

for filename in sys.argv[2:]:
    outfp = open(filename[:-3] + '.graphsize', 'w')
    
    n_kept = 0

    for n, record in enumerate(screed.fastq.fastq_iter(gzip.open(filename))):
        kmer = record['sequence'][:K]
        size = ht.calc_connected_graph_size(kmer, THRESHOLD)
        if size >= THRESHOLD:
            print >>outfp, ">%s\n%s" % (record['name'], record['sequence'])
            n_kept += 1

        if n % 100000 == 0:
            print '...', n, n_kept, n - n_kept + 1
