import sys
import screed.fasta
import os
import khmer
import threading
import Queue
import gc

K = 31                                  # use K-1 for assembly K
HASHTABLE_SIZE = int(4e9)
N_HT = 4

###

MAX_DEGREE = 4

###

WORKER_THREADS = 8
GROUPSIZE = 100


class SequenceGroup(object):
    def __init__(self, order, seqlist):
        self.order = order
        self.seqlist = seqlist


def process(inq, outq, ht):
    global worker_count

    while not done or not inq.empty():
        try:
            g = inq.get(True, 1)
        except Queue.Empty:
            continue

        x = []
        last_record = None
        for record in g.seqlist:
            name = record['name']
            seq = record['sequence']
            trim_seq, trim_at = ht.trim_on_degree(seq, MAX_DEGREE)

            if trim_at > K:
                x.append(record)

        y = [(record['name'], record['sequence']) for record in x]

        gg = SequenceGroup(g.order, y)
        outq.put(gg)

        gc.collect()

    worker_count -= 1


def write(outq, outfp):
    global worker_count
    groups = {}
    next_group = 0
    while worker_count > 0 or not outq.empty():
        try:
            g = outq.get(True, 1)
        except Queue.Empty:
            continue

        groups[g.order] = g

        while next_group in groups:
            g = groups[next_group]
            for name, seq in g.seqlist:
                outfp.write('>%s\n%s\n' % (name, seq,))

            del groups[next_group]
            next_group += 1

            gc.collect()

        if len(groups) > 20:
            print 'WAITAMINIT: len(groups) is', len(groups)


def main():
    global ht, done, worker_count
    done = False
    worker_count = 0

    repfile = sys.argv[1]
    infile = sys.argv[2]
    outprefix = sys.argv[3]

    lowfile = outprefix + '.low'
    highfile = outprefix + '.high'

    print 'saving low-density to:', lowfile
    print 'saving high-density to:', highfile

    print 'making hashtable'
    ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)

    lowfp = open(lowfile, 'w')
    # highfp = open(highfile, 'w')

    print 'eating', repfile
    ht.consume_fasta(repfile)

    inqueue = Queue.Queue(50)
    outqueue = Queue.Queue(50)

    ## worker and writer threads
    for i in range(WORKER_THREADS):
        t = threading.Thread(target=process, args=(inqueue, outqueue, ht))
        worker_count += 1
        t.start()

    threading.Thread(target=write, args=(outqueue, lowfp)).start()

    ### main thread
    x = []
    i = 0
    group_n = 0
    for n, record in enumerate(screed.fasta.fasta_iter(open(infile),
                                                       parse_description=False)):
        if n % 10000 == 0:
            print '...', n

        i += 1
        if i > GROUPSIZE:
            x.append(record)

            g = SequenceGroup(group_n, x)
            inqueue.put(g)
            x = []

            group_n += 1
            i = 0
        else:
            x.append(record)

    # submit last set of sequences
    g = SequenceGroup(group_n, x)
    inqueue.put(g)

    done = True

main()
