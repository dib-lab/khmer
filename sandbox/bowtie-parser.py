#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import sys

'''parses out bowtie map file'''

map_file = open(sys.argv[1], 'r')
out_file = open(sys.argv[2], 'w')
THRESHOLD = 1000
K = 33

for line in map_file:
    data = line.rstrip().split('\t')
    contig_length = int(data[4].split('_')[3]) + K - 1
    read, contig = data[0], data[4]
    if contig_length >= 1000:
        print >>out_file, line.rstrip()
