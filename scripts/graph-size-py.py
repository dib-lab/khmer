import khmer
import sys
import screed
import glob
import os
import subprocess
import zlib
import gzip

K = 32
HASHTABLE_SIZE=4**17+1
THRESHOLD=100

print 'creating ht with size ' + str(HASHTABLE_SIZE)
ht = khmer.new_hashbits(K, HASHTABLE_SIZE, 1)

read_count = 0

fd = open(sys.argv[1], 'w')

files = glob.glob(sys.argv[2])

for filename in files:
    outfp = open(filename + '.graphsize', 'w')
    
    print 'processing file: ' + filename + ' reads processed: ' + str(read_count)
    n_kept = 0
    for n, record in enumerate(screed.fastq.fastq_iter(gzip.open(filename))):
        read_count += 1
        ht.consume(record['sequence'])

        if read_count % 10000 == 0:
            fd.write(str(read_count) + " " + str(ht.n_occupied()) + " " + str(ht.n_occupied() / float(HASHTABLE_SIZE)) +'\n')
            fd.flush()

        kmer = record['sequence'][:K]
        size = ht.calc_connected_graph_size(kmer, THRESHOLD)
        if size >= THRESHOLD:
            print >>outfp, ">%s\n%s" % (record['name'], record['sequence'])
            n_kept += 1

        if n % 10000 == 0:
            print '...', n, n_kept, n - n_kept + 1
