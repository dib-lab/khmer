#! /usr/bin/env python
#
# This script is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import sys

total = 0
totaltotal = int(sys.argv[2])
for line in open(sys.argv[1]):
    size, num, sumnum = line.split()
    size = int(size)
    num = int(num)
    total += num*size

    print size, num, sumnum, total, float(total) / float(totaltotal) * 100.
