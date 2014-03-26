#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import sys
import math
import matplotlib
from scipy.special import gammaln
import pylab

def bloom_mem(n, p):
   n = float(n)
   p = float(p)
   return math.ceil(0 - (n * math.log(p)) / (math.log(2)**2))

def poss_kmers(k):
   '''
   WARNING: This only works with odd numbered kmers due to self-reverse 
   complements in even-numbered kmers (e.g. ACGT).
   '''
   return (4**k) / 2

def inf_mem(k_size, no_kmers):
   k = no_kmers
   n = poss_kmers(k_size)
   x = gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1.0)
   return x * math.log(math.e, 2)

xs = [1000000,
      10000000,
      100000000,
      1000000000,
      10000000000]
bl1 = []
bl5 = []
bl10 = []
bl15 = []
ex19 = []
ex23 = []
ex27 = []
ex31 = []

for x in xs:
   bl_1 = bloom_mem(x, 0.01) / float(x)
   bl_5 = bloom_mem(x, 0.05) / float(x)
   bl_10 = bloom_mem(x, 0.10) / float(x)
   bl_15 = bloom_mem(x, 0.15) / float(x)
   ex_19 = inf_mem(19, x) / float(x)
   ex_23 = inf_mem(23, x) / float(x)
   ex_27 = inf_mem(27, x) / float(x)
   ex_31 = inf_mem(31, x) / float(x)

   bl1.append(bl_1)
   bl5.append(bl_5)
   bl10.append(bl_10)
   bl15.append(bl_15)
   ex19.append(ex_19)
   ex23.append(ex_23)
   ex27.append(ex_27)
   ex31.append(ex_31)

p1, = pylab.semilogx(xs, bl1, 'r--')
p2, = pylab.semilogx(xs, bl5, 'b--')
p3, = pylab.semilogx(xs, bl10, 'g--')
p4, = pylab.semilogx(xs, bl15, 'k--')
p5, = pylab.semilogx(xs, ex19, 'c')
p6, = pylab.semilogx(xs, ex23, 'm')
p7, = pylab.semilogx(xs, ex27, 'y')
p8, = pylab.semilogx(xs, ex31, 'k')

legend_labels = ['Bloom 1%',
                 'Bloom 5%', 
                 'Bloom 10%', 
                 'Bloom 15%',
                 'Exact K=19',
                 'Exact K=23',
                 'Exact K=27',
                 'Exact K=31']

pylab.legend([p1, p2, p3, p4, p5, p6, p7, p8], legend_labels) 
pylab.xlabel("Number of Kmers (Log)")
pylab.ylabel("Bits Needed Per K-mer")

ax = pylab.gca()

for tick in ax.xaxis.get_major_ticks():
   tick.label1.set_fontsize(16)
for tick in ax.yaxis.get_major_ticks():
   tick.label1.set_fontsize(16)

pylab.savefig('memusg.pdf')
