#! /usr/bin/env python
import khmer

filename = 'foo.fa'
fp = open('stats.txt', 'w')

primes = [50000017,
          50000021,
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
          50000347,
          50000387,
          50000389,
          50000393,
          50000399,
          50000429,
          50000449,
          50000477,
          50000491,
          50000507,
          50000513,
          50000537,
          50000539,
          50000557,
          50000563,
          50000567,
          50000569,
          50000617,
          50000627,
          50000641,
          50000683,
          50000711,
          50000719]

_primes = [100000007,
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

_primes = [75000007,
          75000017,
          75000031,
          75000047,
          75000071,
          75000083,
          75000097,
          75000103,
          75000113,
          75000143,
          75000157,
          75000169,
          75000173,
          75000179,
          75000181,
          75000187,
          75000197,
          75000227,
          75000241]

primes = [60000011,
          60000013,
          60000023,
          60000047,
          60000049,
          60000067,
          60000071,
          60000091,
          60000103,
          60000113,
          60000131,
          60000209,
          60000211,
          60000229,
          60000251,
          60000253,
          60000257,
          60000277,
          60000313,
          60000337,
          60000349]

primes = [ 65000011,
           65000021,
           65000029,
           65000051,
           65000057,
           65000063,
           65000071,
           65000081,
           65000119,
           65000137,
           65000141,
           65000147,
           65000161,
           65000167,
           65000179,
           65000209,
           65000213,
           65000233,
           65000249,
           65000261,
           65000269,
           65000303,
           65000321,
           65000357,
           65000393,
           65000399,
           65000437,
           65000443,
           65000473,
           65000549,
           65000561,
           65000567,
           65000569,
           65000597,
           65000603,
           65000623,
           65000633]

###

this_filename = filename
readmask = None

last_n_kept = None
for n, prime in enumerate(primes):
    ht = khmer.new_hashtable(15, prime)

    if not readmask:
        x = ht.consume_fasta_build_readmask(this_filename)
        total_reads, n_consumed, readmask = x

    else:
        total_reads, n_consumed = ht.consume_fasta(this_filename, 0, 0,
                                                   readmask)
        
    print '%d: ate %d k-mers of %d reads' % (n, n_consumed, total_reads)

    print 'filtering...'
    minmax = ht.fasta_file_to_minmax(this_filename, total_reads)
    readmask = ht.filter_fasta_file_any(this_filename, minmax,
                                        5, readmask)

    n_seq_kept = readmask.n_kept()

    print '%d: kept %d of %d (%.1f%%)' % (n, n_seq_kept, total_reads, n_seq_kept/float(total_reads)*100)

    fp.write('%d %d %f\n' % (n, n_seq_kept, ht.n_occupied() / float(prime) ))
    fp.flush()

    if n_seq_kept == last_n_kept:
        print 'DONE'
        break

    last_n_kept = n_seq_kept

###

