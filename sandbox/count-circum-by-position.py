import sys, screed.fasta, os
import khmer

K = 32
HASHTABLE_SIZE=int(8e9)
N_HT = 4

###

RADIUS=2
MAX_CIRCUM=4                            # 4 seems to eliminate lump in 1m.fa
MAX_VOLUME=200

infile = sys.argv[1]
outfile = sys.argv[2]
outfp = open(outfile, 'w')
    
print 'saving to:', outfile

print 'making hashtable'
ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)

print 'eating', infile
ht.consume_fasta(infile)

hist = [0.0] * 200
histcount = [0] * 200

for n, record in enumerate(screed.fasta.fasta_iter(open(infile))):
    if n % 10000 == 0:
        print '... saving', n
    seq = record['sequence']

    for pos in range(0, len(seq) - K):
        circum = ht.count_kmers_within_radius(seq[pos:pos+K], RADIUS,
                                              MAX_VOLUME)
        hist[pos] += circum
        histcount[pos] += 1

for i in range(len(hist)):
    if histcount[i]:
        print >>outfp, i, hist[i], histcount[i], hist[i]/histcount[i]
