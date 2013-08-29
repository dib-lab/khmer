import sys
import screed.fasta
import os
import khmer
from khmer.thread_utils import ThreadedSequenceProcessor, verbose_fastq_iter

K = 32
HT_SIZE = 4e9
N_HT = 4

WORKER_THREADS = 8
GROUPSIZE = 100

###


def main():
    repfile = sys.argv[1]
    infile = sys.argv[2]

    outfile = os.path.basename(infile) + '.fno255'
    if len(sys.argv) >= 4:
        outfile = sys.argv[3]

    print 'file to count from: %s' % repfile
    print 'input file to filter: %s' % infile
    print 'filtering to output:', outfile
    print '-- settings:'
    print 'K', K
    print 'N THREADS', WORKER_THREADS
    print '--'

    print 'making hashtable'
    ht = khmer.new_counting_hash(K, HT_SIZE, N_HT)

    print 'consuming input', repfile
    ht.consume_fasta(repfile)

    outfp = open(outfile, 'w')

    def process_fn(record, ht=ht):
        name = record['name']
        seq = record['sequence']
        if 'N' in seq:
            return None, None

        if len(seq) < K:
            return None, None

        if ht.get_max_count(seq) >= 255:
            return None, None

        return name, seq

    tsp = ThreadedSequenceProcessor(process_fn, WORKER_THREADS, GROUPSIZE)

    ###

    tsp.start(verbose_fastq_iter(infile), outfp)

if __name__ == '__main__':
    main()
