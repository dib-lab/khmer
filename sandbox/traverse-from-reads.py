import khmer
import sys
import threading
import Queue
import gc
import os.path

K = 20
HASHTABLE_SIZE = int(1e9)
N_HT = 4

COUNTING_SIZE = int(1e8)

ht = None

###


def main(filename):
    global ht

    basename = os.path.basename(filename)

    print 'input file to partition: %s' % filename
    print '-- settings:'
    print 'K', K
    print 'HASHTABLE SIZE %g' % HASHTABLE_SIZE
    print 'N HASHTABLES %d' % N_HT
    print '--'

    ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)

    ht.consume_fasta(filename)

    counting = khmer.new_counting_hash(K, COUNTING_SIZE, N_HT)
    ht.traverse_from_reads(filename, 100, 5000, 5, counting)

    print 'saving stoptags binary'
    ht.save_stop_tags(basename + '.stoptags')
    print 'saving stoptags text'
    ht.print_stop_tags(basename + '.stoptags.txt')

    sys.exit(0)

if __name__ == '__main__':
    main(sys.argv[1])
