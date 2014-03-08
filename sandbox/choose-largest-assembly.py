#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
import sys
from screed.fasta import fasta_iter


def count_sum_contigs(cutoff, filename):
    total = 0
    for record in fasta_iter(open(filename)):
        seqlen = len(record['sequence'])
        if seqlen >= cutoff:
            total += seqlen

    return total

min_length = int(sys.argv[1])

best_filename = sys.argv[2]
best_sum = count_sum_contigs(min_length, best_filename)
print >>sys.stderr, 'at:', best_sum, best_filename

for filename in sys.argv[3:]:
    this_sum = count_sum_contigs(min_length, filename)

    print >>sys.stderr, 'at:', this_sum, filename
    if this_sum > best_sum:
        best_sum = this_sum
        best_filename = filename

if best_sum == 0:
    print >>sys.stderr, "no non-zero assembly"
    sys.exit(0)

print >>sys.stderr, 'keeping:', best_sum, best_filename

for record in fasta_iter(open(best_filename)):
    if len(record['sequence']) >= min_length:
        print '>%s\n%s' % (record['name'], record['sequence'],)
