#! /usr/bin/env python
import khmer, sys

def report_and_checkpoint(name, count1, count2):
    global ht
    print name, count1, count2

    if name == 'do_truncated_partition/read' and count1 % 100000 == 0:
        ht.save_checkpoint(sys.argv[1] + '.pmap',
                           sys.argv[1] + '.surrender')

K=32
HASHTABLE_SIZE=int(4**16)+1

infile = sys.argv[1]
outfile = sys.argv[2]

ht = khmer.new_hashtable(K, HASHTABLE_SIZE)

#ht.load_checkpoint('x/1m-filtered.fa.pmap.end', 'x/1m-filtered.fa.surrender.end')

#n_partitions = ht.do_exact_partition(infile)
n_partitions = ht.do_truncated_partition(infile, outfile, report_and_checkpoint)
print n_partitions, 'partitions kept'

ht.save_checkpoint(sys.argv[1] + '.pmap.end',
                   sys.argv[1] + '.surrender.end')
