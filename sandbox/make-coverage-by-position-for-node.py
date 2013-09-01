#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import sys
import csv
import matplotlib
from pylab import *
import numpy

filein = open(sys.argv[1], 'r')
contig_seq = open(sys.argv[2], 'r')

reader = csv.reader(filein)


xlist = []
ylist = []
stat_list = []
for n, row in enumerate(reader):
    xlist.append(n)
    ylist.append(int(row[0]))
    if int(row[0]) > 0:
        stat_list.append(int(row[0]))

average = numpy.average(stat_list)
sd = numpy.std(stat_list)


n2_list = []
scaffold_list = []
for n, line in enumerate(contig_seq):
    if n == 1:
        line = list(line.rstrip())
        for n2, each in enumerate(line):
            if each == 'N':
                n2_list.append(int(n2) - 32 + 1)
                scaffold_list.append(1)
            else:
                n2_list.append(int(n2) - 32 + 1)
                scaffold_list.append(0)

figure()
plot(xlist, ylist, 'r-')
plot(n2_list, scaffold_list, 'b-')
titletext = '%s avg %.2f sd %.2f' % (sys.argv[1], average, sd)
title(titletext)
# show()
savefig(sys.argv[1] + '.png')

print 'average=', numpy.average(stat_list)
print 'std dev=', numpy.std(stat_list)
