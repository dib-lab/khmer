#! /usr/bin/env python
#
# This script is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import sys

'''prints out summary of all stats for all partition stats in a directory'''

fp = open(sys.argv[1], 'w')
stat_files = sys.argv[2:]

for file in stat_files:
    f = open(file, 'r')
    for n, line in enumerate(f):
        if n == 2:
            line = line.rstrip().split()
            total_contigs = line[0]
            length = line[1]
            max_contig = line[2]
            group_index = line[3].find('group')
            group = line[3][group_index + 5:group_index + 9]
            print >>fp, '%s %s %s %s' % (
                group, total_contigs, length, max_contig)
    f.close()
