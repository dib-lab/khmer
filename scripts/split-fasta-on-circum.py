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
outprefix = sys.argv[2]

lowfile = outprefix + '.low'
highfile = outprefix + '.high'

print 'saving low-density to:', lowfile
print 'saving high-density to:', highfile

print 'making hashtable'
ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)

lowfp = open(lowfile, 'w')
highfp = open(highfile, 'w')

print 'eating', infile
ht.consume_fasta(infile)

for n, record in enumerate(screed.fasta.fasta_iter(open(infile),
                                                   parse_description=False)):
    if n % 10000 == 0:
        print '... saving', n

    seq = record['sequence']

    # calculate circumference for every point.
    end = len(seq) - K
    is_high = False

    for pos in range(end, -1, -1):
        circum = ht.count_kmers_on_radius(seq[pos:pos+K], RADIUS, MAX_VOLUME)
        print seq[pos:pos+K], circum

        if circum >= MAX_CIRCUM:
            is_high = True
            break

    # sort "high circumference" and "low" circumerence sequences separately.
    if is_high:
        fp = highfp
    else:
        fp = lowfp
        
    print >>fp, '>%s\n%s' % (record['name'], seq)
