import khmer
import sys
import screed
import os.path
import threading, Queue

K = 32
HASHTABLE_SIZE=int(4e9)
THRESHOLD=100
N_HT=4

###

RADIUS=2
MAX_CIRCUM=4                            # 4 seems to eliminate lump in 1m.fa
MAX_VOLUME=200

incr = 2*RADIUS

###

GROUPSIZE=5
WORKER_THREADS=4

###

def is_pair(r1, r2):
    a = r1['name'].split('/')[0]
    b = r2['name'].split('/')[0]

    return (a==b)

def trim_by_circumference(ht, name, seq):
    # calculate circumference for every point.
    end = len(seq) - K
    is_high = False

    pos = 0
    for pos in range(0, end, incr):
        circum = ht.count_kmers_on_radius(seq[pos:pos+K], RADIUS, MAX_VOLUME)

        if circum >= MAX_CIRCUM:
            is_high = True
            break

    # ok. sequence has high-radius k-mers; can we trim them off?
    if is_high and pos > incr:
        pos -= incr

        # find last k-mer with a low radius:
        i = 1
        for i in range(1, incr):
            circum = ht.count_kmers_on_radius(seq[pos+i:pos+i+K],
                                              RADIUS, MAX_VOLUME)
            if circum >= MAX_CIRCUM:
                break

        pos += i - 1

        # now trim sequence:
        seq = seq[:pos+K]
        is_high = False
        name += "\tTRUNC.%d" % pos

    if is_high:
        return (None, None)
    
    return (name, seq)

def process(inq, outq, ht):
    global worker_count
    
    while not done or not inq.empty():
        try:
            recordlist = inq.get(True, 1)
        except Queue.Empty:
            continue

        x = []
        last_record = None
        for record in recordlist:
            kmer = record['sequence'][:K]
            size = ht.calc_connected_graph_size(kmer, THRESHOLD, True)
            if size >= THRESHOLD:
                # keep pairs together if either is "good"
                if last_record and is_pair(last_record, record):
                    x.append(last_record)
                x.append(record)
                record = None

            last_record = record

        y = []
        for record in x:
            name, seq = trim_by_circumference(ht, record['name'],
                                              record['sequence'])
            if name:
                y.append((name, seq))

        outq.put(y)

    worker_count -= 1

def write(outq, outfp):
    global worker_count
    while worker_count > 0 or not outq.empty():
        try:
            recordlist = outq.get(True, 1)
        except Queue.Empty:
            continue

        for name, seq in recordlist:
            outfp.write('>%s\n%s\n' % (name, seq,))

def main():
    global done, worker_count
    done = False
    worker_count = 0
    
    infile = sys.argv[1]
    outfile = os.path.basename(infile) + '.graphcirc'
    if len(sys.argv) == 3:
        outfile = sys.argv[2]

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
    for n, record in enumerate(screed.fasta.fasta_iter(open(infile))):
        if n % 10000 == 0:
            print '...', n

        i += 1
        if i > GROUPSIZE:
            this_name = record['name'].split('/')[0]
            last_name = x[-1]['name'].split('/')[0]

            if is_pair(record, x[-1]):     # preserve pairs
                x.append(record)
                inqueue.put(x)
                x = []
            else:
                inqueue.put(x)
                x = [record]
                
            i = 0
        else:
            x.append(record)

    # submit last set of sequences
    inqueue.put(x)

    done = True

if __name__ == '__main__':
    main()
