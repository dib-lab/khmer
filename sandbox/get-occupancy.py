#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import khmer
import sys
import screed
import os
import subprocess
import glob
import gzip

K = 32
HASHTABLE_SIZE = int(1e9)
# HASHTABLE_SIZE = 1000000000

print 'creating ht with size ' + str(HASHTABLE_SIZE)
ht = khmer.new_hashbits(K, HASHTABLE_SIZE, 1)

read_count = 0

files = sys.argv[1:]

for filename in files:
    print "processing file: " + filename + " reads processed: " + \
        str(read_count)

    for n, record in enumerate(screed.fasta.fasta_iter(open(filename))):
        read_count += 1
        if len(record['sequence']) >= K:
            ht.consume(record['sequence'])

        if read_count % 10000 == 0:
            print str(read_count), str(ht.n_occupied())

print str(read_count), str(ht.n_occupied())

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
