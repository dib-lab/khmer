import sys
import screed.fasta
import os
import khmer

K = 32
HASHTABLE_SIZE = int(8e9)
N_HT = 4

MAX_RADIUS = 500
VOLUME = 500

###

infile = sys.argv[1]
outfile = sys.argv[2]
if len(sys.argv) > 3:
    RADIUS = int(sys.argv[3])

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

    middle = (len(seq) - K + 1) / 2

    r = ht.find_radius_for_volume(seq[middle:middle + K], VOLUME, MAX_RADIUS)

    print >>outfp, '>%s r=%d\n%s' % (record['name'], r, record['sequence'])
