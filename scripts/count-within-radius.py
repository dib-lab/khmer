import sys, screed.fasta, os
import khmer

K = 32
HASHTABLE_SIZE=int(8e9)
N_HT = 4
THRESHOLD=100

infile = sys.argv[1]
outfile = sys.argv[2]
print 'saving to:', outfile

print 'making hashtable'
ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)

print 'eating', infile
ht.consume_fasta(infile)

outfp = open(outfile, 'w')
for n, record in enumerate(screed.fasta.fasta_iter(open(infile))):
    if n % 10000 == 0:
        print '... saving', n
    seq = record['sequence']
    
    area = ht.count_kmers_within_radius(seq[:K], THRESHOLD, 2000)
    print >>outfp, '>%s a=%d\n%s' % (record['name'], area, record['sequence'])
