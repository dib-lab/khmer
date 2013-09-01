#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import sys
import math

#
# the two unknowns that we know for simulations but want to calculate for
# unknown situations.
#

P_ERROR = .01                           # per-base error rate (an unknown)

t=float(99366)                          # here, T is the total # of unique
                                        # k-mers (an unknown, normally)

# read in the data from 'swr.py'.
lines = open(sys.argv[1]).readlines()
lines = [ x.split() for x in lines ]
lines = [ (int(x[0]), int(x[1])) for x in lines ]

# here, 'x' is the number of k-mers sampled, 'y' is the number of unique
# k-mers seen:
xdata = [ x[0] for x in lines ]
ydata = [ x[1] for x in lines ]

# calculate the expected curve.
for x in xdata:
    y2 = T*(1. - math.exp(-x / T)) + x*P_ERROR*3./4.
    print x, y2
