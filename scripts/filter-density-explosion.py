import sys, screed.fasta, os
import khmer
from khmer.thread_utils import ThreadedSequenceProcessor, verbose_fasta_iter

K = 32
HASHTABLE_SIZE=int(8e8)
N_HT = 4

RADIUS=2
MAX_VOLUME=10

WORKER_THREADS=8
GROUPSIZE=100

###

def main():
    repfile = sys.argv[1]
    infile = sys.argv[1]
    if len(sys.argv) >= 3:
        infile = sys.argv[2]
        
    outfile = os.path.basename(infile) + '.loess'
    if len(sys.argv) >= 4:
        outfile = sys.argv[3]

    print 'file with representative artifacts: %s' % repfile
    print 'input file to degree filter: %s' % infile
    print 'filtering to output:', outfile
    print '-- settings:'
    print 'K', K
    print 'HASHTABLE SIZE %g' % HASHTABLE_SIZE
    print 'N HASHTABLES %d' % N_HT
    print 'N THREADS', WORKER_THREADS
    print 'RADIUS', RADIUS
    print 'MAX DENSITY', MAX_VOLUME / RADIUS
    print '--'

    print 'making hashtable'
    ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)

    outfp = open(outfile, 'w')

    print 'eating', repfile
    ht.consume_fasta(repfile)

    def process_fn(record, ht=ht):
        name = record['name']
        seq = record['sequence']
        if 'N' in seq:
            return None, None

        trim_seq, trim_at = ht.trim_on_density_explosion(seq, RADIUS,
                                                         MAX_VOLUME)

        if trim_at >= K:
            return name, trim_seq

        return None, None

    tsp = ThreadedSequenceProcessor(process_fn, WORKER_THREADS, GROUPSIZE)

    ###

    tsp.start(verbose_fasta_iter(infile), outfp)

if __name__ == '__main__':
    main()
