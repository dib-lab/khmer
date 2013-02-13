#! /usr/bin/env python
import khmer

filename = '1m.fa'

ht_full = khmer.new_hashtable(15, 4 ** 15)
print 'filtering 15-mers exactly:'
total_reads, n_consumed = ht_full.consume_fasta(filename)
print 'ate %d k-mers of %d reads' % (n_consumed, total_reads)
print 'filtering...'
if 0:
    (total_reads_2, n_seq_kept) = khmer.filter_fasta_file(ht_full,
                                                          filename,
                                                          total_reads,
                                                          filename +
                                                          '.15.exact',
                                                          5)

    print 'kept %d of %d (%.1f%%)' % (
        n_seq_kept, total_reads, n_seq_kept / float(total_reads) * 100)

print 'counting!'
print '%d total k-mers' % (ht_full.n_occupied())

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
