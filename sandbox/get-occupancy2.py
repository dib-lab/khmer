#! /usr/bin/env python
#
# This script is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import khmer
import sys
import os.path

K = 32
HASHTABLE_SIZE = int(1e9)

ht = khmer.new_hashbits(K, HASHTABLE_SIZE, 1)


def callback(name, total_reads, n_consumed):
    print name, total_reads, n_consumed, ht.n_occupied()


def main(filename):
    basename = os.path.basename(filename)

    # populate the hash table and tag set
    ht.consume_fasta(filename, 0, 0, None, True, callback)

if __name__ == '__main__':
    main(sys.argv[1])
