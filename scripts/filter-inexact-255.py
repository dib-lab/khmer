#! /usr/bin/env python
import sys

import khmer

KSIZE=20
HASHTABLE_SIZE=int(2e9)
MIN_ABUNDANCE=255

infilename = sys.argv[1]
outfilename = sys.argv[2]

print 'filtering %d-mers exactly;' % KSIZE
print 'from file %s to file %s;' % (infilename, outfilename,)
print 'looking for sequences containing k-mers with all counts < %d' \
      % MIN_ABUNDANCE

ht = khmer.new_counting_hash(KSIZE, HASHTABLE_SIZE, 4)
total_reads, n_consumed = ht.consume_fasta(infilename)

print 'just ate %d reads, %d k-mers' % (total_reads, n_consumed)

print 'filtering...'

minmax = ht.fasta_file_to_minmax(infilename, total_reads)
readmask = ht.filter_fasta_file_any(minmax, MIN_ABUNDANCE)
readmask.invert()

print 'keeping %d reads' % readmask.n_kept()

print 'writing...'
readmask.filter_fasta_file(infilename, outfilename)

print 'hashtable occupancy:', ht.n_occupied() / float(HASHTABLE_SIZE)
