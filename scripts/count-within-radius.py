import sys, screed.fasta, os
import khmer

K = 32
HASHTABLE_SIZE=int(8e9)
N_HT = 4
RADIUS=10

###

MAX_DENSITY=2000

infile = sys.argv[1]
outfile = sys.argv[2]
if len(sys.argv) > 3:
    RADIUS=int(sys.argv[3])
    
print 'saving to:', outfile

print 'making hashtable'
ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)

print 'eating', infile
ht.consume_fasta(infile)

outfp = open(outfile, 'w')
for n, record in enumerate(screed.open(infile)):
    if n % 10000 == 0:
        print '... saving', n
    seq = record['sequence']

    middle = (len(seq) - K + 1) / 2
    
    density = ht.count_kmers_within_radius(seq[middle:middle+K], RADIUS,
                                           MAX_DENSITY)
    density /= float(RADIUS)
    
    print >>outfp, '>%s d=%.3f\n%s' % (record['name'], density, record['sequence'])
