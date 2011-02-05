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
circfile = outprefix + '.circ'

print 'saving low-density to:', lowfile
print 'saving high-density to:', highfile
print 'saving all densities to:', circfile

print 'making hashtable'
ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)

lowfp = open(lowfile, 'w')
highfp = open(highfile, 'w')
circfp = open(circfile, 'w')

print 'eating', infile
ht.consume_fasta(infile)

start = RADIUS
incr = 2*RADIUS

for n, record in enumerate(screed.fasta.fasta_iter(open(infile),
                                                   parse_description=False)):
    if n % 10000 == 0:
        print '... saving', n

    seq = record['sequence']
    end = len(seq) - K + 1 - incr/2

    is_high = False
    circums = []
    for pos in range(end, start, -incr):
        circum = ht.count_kmers_on_radius(seq[pos:pos+K], RADIUS, MAX_VOLUME)

        circums.append((circum, pos))
        print >>circfp, circum
        
        if circum >= MAX_CIRCUM:
            is_high = True

    if is_high:
        for (circum, pos) in circums:
            if circum < MAX_CIRCUM:
                break

        chop_end = pos
            
        circums.reverse()
        for (circum, pos) in circums:
            if circum < MAX_CIRCUM:
                break
        chop_start = pos

        if max([ circ for (circ,_) in circums ]) < MAX_CIRCUM and \
               chop_end - chop_start >= 0:
            sequence = record['sequence'][chop_start:chop_end + K]
            record['name'] = record['name'] + '\tTRUNC:%d-%d' % (chop_start,
                                                                 chop_end + K)
            is_high = False

    if is_high:
        fp = highfp
    else:
        fp = lowfp
        
    print >>fp, '>%s\n%s' % (record['name'], record['sequence'])
