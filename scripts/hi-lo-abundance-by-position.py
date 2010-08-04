import sys
import os
import khmer

HASHTABLE_SIZE=4**12+1
K=32

def write_dist(dist, fp):
    for n, i in enumerate(dist):
        fp.write('%d %d\n' % (n, i))

filename = sys.argv[1]
outfile = os.path.basename(filename)

ht = khmer.new_hashtable(K, HASHTABLE_SIZE)
ht.consume_fasta(filename)
x = ht.fasta_count_kmers_by_position(filename, 100, 1)
write_dist(x, open(outfile + '.pos.abund=1', 'w'))

y = ht.fasta_count_kmers_by_position(filename, 100, 255)
write_dist(y, open(outfile + '.pos.abund=255', 'w'))
