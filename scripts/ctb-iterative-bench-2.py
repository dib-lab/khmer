#! /usr/bin/env python
import khmer

filename = '1m.fa'
fp = open('stats.txt', 'w')

primes = [50000021,
          50000047,
          50000059,
          50000063,
          50000101,
          50000131,
          50000141,
          50000161,
          50000201,
          50000207,
          50000221,
          50000231,
          50000233,
          50000239,
          50000243,
          50000257,
          50000273,
          50000309,
          50000311,
          50000329,
          50000341,
          50000347]

primes = [100000007,
          100000037,
          100000039,
          100000049,
          100000073,
          100000081,
          100000123,
          100000127,
          100000193,
          100000213,
          100000217,
          100000223,
          100000231,
          100000237,
          100000259,
          100000267,
          100000279,
          100000357,
          100000379]

###

this_filename = filename
for n, prime in enumerate(primes):
    ht_small = khmer.new_hashtable(15, prime)
    next_filename = filename + '.15.round%d' % n

    total_reads, n_consumed = ht_small.consume_fasta(this_filename)
    print '%d: ate %d k-mers of %d reads' % (n, n_consumed, total_reads)

    print 'filtering...'
    (total_reads_2, n_seq_kept) = khmer.filter_fasta_file(ht_small,
                                                          this_filename,
                                                          total_reads,
                                                          next_filename,
                                                          5)

    print '%d: kept %d of %d (%.1f%%)' % (n, n_seq_kept, total_reads, n_seq_kept/float(total_reads)*100)

    fp.write('%d %d\n' % (n, n_seq_kept))
    fp.flush()

    this_filename = next_filename

###

