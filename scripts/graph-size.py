import khmer
import sys
import screed
import os.path
import threading, Queue

K = 32
HASHTABLE_SIZE=int(4e9)
THRESHOLD=500
N_HT=4
WORKER_THREADS=5

###

GROUPSIZE=100

###

class SequenceGroup(object):
    def __init__(self, order, seqlist):
        self.order = order
        self.seqlist = seqlist

def is_pair(r1, r2):
    a = r1['name'].split('/')[0]
    b = r2['name'].split('/')[0]

    return (a==b)

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
            kmer = record['sequence'][:K]
            size = ht.calc_connected_graph_size(kmer, THRESHOLD)
            if size >= THRESHOLD:
                # keep pairs together if either is "good"
                if last_record and is_pair(last_record, record):
                    x.append(last_record)
                x.append(record)
                record = None

            last_record = record

        y = [ (record['name'], record['sequence']) for record in x ]

        gg = SequenceGroup(g.order, y)
        outq.put(gg)

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

def main():
    global done, worker_count
    done = False
    worker_count = 0
    
    infile = sys.argv[1]
    outfile = os.path.basename(infile) + '.graphsize'
    if len(sys.argv) == 3:
        outfile = sys.argv[2]

    print 'input file to graphsize filter: %s' % infile
    print 'filtering to output:', outfile
    print '-- settings:'
    print 'K', K
    print 'HASHTABLE SIZE %g' % HASHTABLE_SIZE
    print 'N HASHTABLES %d' % N_HT
    print 'THRESHOLD', THRESHOLD
    print 'N THREADS', WORKER_THREADS
    print '--'

    print 'creating ht'
    ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)
    print 'eating fa', infile
    total_reads, n_consumed = ht.consume_fasta(infile)
    outfp = open(outfile, 'w')

    inqueue = Queue.Queue(50)
    outqueue = Queue.Queue(50)

    ## worker and writer threads
    for i in range(WORKER_THREADS):
        t = threading.Thread(target=process, args=(inqueue, outqueue, ht))
        worker_count += 1
        t.start()

    threading.Thread(target=write, args=(outqueue, outfp)).start()

    ### main thread
    x = []
    i = 0
    group_n = 0
    for n, record in enumerate(screed.fasta.fasta_iter(open(infile))):
        if n % 10000 == 0:
            print '...', n

        i += 1
        if i > GROUPSIZE:
            this_name = record['name'].split('/')[0]
            last_name = x[-1]['name'].split('/')[0]

            if is_pair(record, x[-1]):     # preserve pairs
                x.append(record)

                g = SequenceGroup(group_n, x)
                inqueue.put(g)
                x = []
            else:
                g = SequenceGroup(group_n, x)
                inqueue.put(g)
                x = [record]

            group_n += 1
            i = 0
        else:
            x.append(record)

    # submit last set of sequences
    g = SequenceGroup(group_n, x)
    inqueue.put(g)

    done = True

if __name__ == '__main__':
    main()
