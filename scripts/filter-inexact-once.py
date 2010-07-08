#! /usr/bin/env python
import sys

import khmer

KSIZE=32
HASHTABLE_SIZE=4**15
MIN_ABUNDANCE=2

infilename = sys.argv[1]
outfilename = sys.argv[2]

print 'filtering %d-mers exactly;' % KSIZE
print 'from file %s to file %s;' % (infilename, outfilename,)
print 'looking for sequences containing k-mers, all with count >= %d' \
      % MIN_ABUNDANCE

ht = khmer.new_hashtable(KSIZE, HASHTABLE_SIZE)
total_reads, n_consumed = ht.consume_fasta(infilename)

print 'just ate %d reads, %d k-mers' % (total_reads, n_consumed)

print 'filtering...'

readmask = ht.filter_fasta_file_run(infilename, total_reads, MIN_ABUNDANCE,
                                    KSIZE)

print 'keeping %d reads' % readmask.n_kept()

print 'writing...'
readmask.filter_fasta_file(infilename, outfilename)

print 'hashtable occupancy:', ht.n_occupied() / float(HASHTABLE_SIZE)
