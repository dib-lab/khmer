import sys
import screed.fasta
import os
import khmer
from khmer.thread_utils import ThreadedSequenceProcessor, verbose_fasta_iter

K = 31                                  # use K-1 for part/assembly with K
HASHTABLE_SIZE = int(4e9)
N_HT = 4
MAX_DEGREE = 3

WORKER_THREADS = 8
GROUPSIZE = 100

###


def main():
    repfile = sys.argv[1]
    infile = sys.argv[1]
    if len(sys.argv) >= 3:
        infile = sys.argv[2]

    outfile = os.path.basename(infile) + '.low'
    if len(sys.argv) >= 4:
        outprefix = sys.argv[3]

    print 'file with representative artifacts: %s' % repfile
    print 'input file to degree filter: %s' % infile
    print 'filtering to output:', outfile
    print '-- settings:'
    print 'K', K
    print 'HASHTABLE SIZE %g' % HASHTABLE_SIZE
    print 'N HASHTABLES %d' % N_HT
    print 'N THREADS', WORKER_THREADS
    print 'MAX DEGREE', MAX_DEGREE
    print '--'

    print 'making hashtable'
    ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)

    outfp = open(outfile, 'w')

    print 'eating', repfile
    ht.consume_fasta(repfile)

    def process_fn(record, ht=ht):
        name = record['name']
        seq = record['sequence']
        trim_seq, trim_at = ht.trim_on_degree(seq, MAX_DEGREE)

        if trim_at > K:
            return name, trim_seq

        return None, None

    tsp = ThreadedSequenceProcessor(process_fn, WORKER_THREADS, GROUPSIZE)

    ###

    tsp.start(verbose_fasta_iter(infile), outfp)

if __name__ == '__main__':
    main()
