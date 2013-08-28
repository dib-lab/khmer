#! /usr/bin/env python
#
# This script is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import sys
import os
import khmer


def write_dist(dist, fp):
    for n, i in enumerate(dist):
        fp.write('%d %d\n' % (n, i))

hashfile = sys.argv[1]
filename = sys.argv[2]
outfile = os.path.basename(filename)

ht = khmer.new_counting_hash(1, 1, 1)
# ht.consume_fasta(filename)
ht.load(hashfile)

x = ht.fasta_count_kmers_by_position(filename, 100, 1)
write_dist(x, open(outfile + '.pos.abund=1', 'w'))

y = ht.fasta_count_kmers_by_position(filename, 100, 255)
write_dist(y, open(outfile + '.pos.abund=255', 'w'))
