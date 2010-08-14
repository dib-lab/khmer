#! /usr/bin/env python
import khmer, sys, screed, threading

K=32
HASHTABLE_SIZE=int(4**15)+1

PERIOD=100000
CHUNK=10000

infile = sys.argv[1]
outfile = sys.argv[2]
min_partition_size = int(sys.argv[3])

ht = khmer.new_hashtable(K, HASHTABLE_SIZE)

l = []
m = []

def find_tags(ht, kmer_l):
    global m
    for kmer in kmer_l:
        ppi = ht.find_all_tags(kmer)
        m.append(ppi)

def do_assign(ht, l):
    global m
    
    tl = []
    for i in range(0, len(l), CHUNK):
        t = threading.Thread(target=find_tags, args=(ht, l[i:i+CHUNK]))
        tl.append(t)
        t.start()

    print 'started', len(tl), len(l), len(l) / CHUNK
    for t in tl:
        t.join()

    for ppi in m:
        if ppi:
            ht.assign_partition_id_th(ppi)
    m = []

for n, record in enumerate(screed.fasta.fasta_iter(open(infile))):
    if n % 10000 == 0:
        print '...', n
    seq = record['sequence']
    try:
        ht.consume(seq)
        l.append(seq[:K])
    except ValueError:
        pass

    if len(l) == PERIOD:
        do_assign(ht, l)
        l = []

if l:
    do_assign(ht, l)

n_kept = ht.output_partitions(infile, outfile)
print n_kept, 'partitions kept'
