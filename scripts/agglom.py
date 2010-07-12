#! /usr/bin/env python
import sys
import khmer

K=32
HTSIZE = 4**17+1
CUTOFF = 4

def consume_fasta_if_intersect(ht, filename, total_reads, orig_readmask):
    #mmt = ht.fasta_file_to_minmax(filename, total_reads, orig_readmask)
    #new_readmask = ht.filter_fasta_file_any(mmt, 2)
    new_readmask = ht.filter_fasta_file_run(filename, total_reads, 1, 5)
    print 'XXX', new_readmask.n_kept()
    new_readmask.save('the_readmask')
    ht = khmer.new_hashtable(K, HTSIZE)
    (t, n) = ht.consume_fasta(filename, 0, 0, new_readmask, False)
    return ht, n

filename = sys.argv[1]
outfile = sys.argv[2]

ht = khmer.new_hashtable(K, HTSIZE)
(total_reads, n_consumed, readmask_orig) = \
              ht.consume_fasta_build_readmask(filename)

print 'ate %d reads' % total_reads

mmt = ht.fasta_file_to_minmax(filename, total_reads, readmask_orig)

readmask_all = ht.filter_fasta_file_all(mmt, CUTOFF)
readmask_any = ht.filter_fasta_file_any(mmt, CUTOFF)

## ok, now: delete hashtable; create new one; populate with all >= CUTOFF.

ht = khmer.new_hashtable(K, HTSIZE)

# go through and "eat" the all >= CUTOFF
print '---'
(_, n_consumed_1) = ht.consume_fasta(filename, 0, 0, readmask_all, False)

# go through again and eat those that intersect with the 'all >= CUTOFF' reads
print '---'
ht, n_consumed_2 = consume_fasta_if_intersect(ht, filename, total_reads, readmask_any)

print '---'
ht, n_consumed_3 = consume_fasta_if_intersect(ht, filename, total_reads, readmask_orig)

print '---'
ht, n_consumed_4 = consume_fasta_if_intersect(ht, filename, total_reads, readmask_orig)

print '---'
ht, n_consumed_5 = consume_fasta_if_intersect(ht, filename, total_reads, readmask_orig)

print '---'
ht, n_consumed_6 = consume_fasta_if_intersect(ht, filename, total_reads, readmask_orig)

print (n_consumed, n_consumed_1, n_consumed_2, n_consumed_3, n_consumed_4, n_consumed_5, n_consumed_6)
