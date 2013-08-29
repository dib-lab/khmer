import sys
import screed.fasta
import os
import khmer

K = 32
HASHTABLE_SIZE = int(2e9)
N_HT = 4
RADIUS = 4

###

MAX_DENSITY = 2000
THRESHOLD = 10

infile = sys.argv[1]
outfile = sys.argv[2]
if len(sys.argv) > 3:
    RADIUS = int(sys.argv[3])

print 'saving to:', outfile

print 'making hashtable'
ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)

print 'eating', infile
ht.consume_fasta(infile)

seen = set()

outfp = open(outfile, 'w')
for n, record in enumerate(screed.open(infile)):
    if n % 10000 == 0:
        print '... saving', n
        if n > 100000:
            break

    seq = record['sequence']

    for pos in range(0, len(seq) - K + 1):
        kmer = seq[pos:pos + K]
        if kmer in seen:
            continue

        density = ht.count_kmers_within_radius(kmer, RADIUS, MAX_DENSITY)
        if density >= THRESHOLD:
            seen.add(kmer)
            print >>outfp, kmer, density
