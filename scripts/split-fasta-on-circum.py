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

start = RADIUS
incr = 2*RADIUS

for n, record in enumerate(screed.fasta.fasta_iter(open(infile),
                                                   parse_description=False)):
    if n % 10000 == 0:
        print '... saving', n

    seq = record['sequence']

    # between [RADIUS:-RADIUS] kmers, calculate circumference every 2*RADIUS.
    # chop sequence from first high circumference, onwards.
    
    end = len(seq) - K + 1 - incr/2
    is_high = False

    for pos in range(start, end, incr):
        circum = ht.count_kmers_on_radius(seq[pos:pos+K], RADIUS, MAX_VOLUME)

        if circum >= MAX_CIRCUM:
            is_high = True
#            if pos == start:    # entire sequence is crud => high file
#                is_high = True
#            else:
#                chop = pos - incr
#                seq = seq[:chop + K]    # may be salvageable
                
            break

    # sort "high circumference" and "low" circumerence sequences separately.
    if is_high:
        fp = highfp
    else:
        fp = lowfp
        
    print >>fp, '>%s\n%s' % (record['name'], seq)
