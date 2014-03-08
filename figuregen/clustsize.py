#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
import khmer
import sys
import random
import math
import numpy

def callback(a, b, c):
   pass

bases = ['A', 'C', 'G', 'T']

#p = float(sys.argv[1])

def estimate_mean(x):
   ''' 
   uses simple percentile bootstrap to generate confidence intervals
   '''
   n = len(x)
   upper = 0
   med = 0
   lower = 0

   means = []

   for i in range(200):
      tmp = []
      for j in range(n):
         tmp.append(x[random.randint(0,n-1)])

      means.append(numpy.mean(tmp))

   means = sorted(means)

   lower = means[4]
   med = means[99]
   upper = means[194]

   return lower, med, upper

def calc_ht_size(m, k):
   return int(m / k)

def calc_m(n, p):
   n = float(n)
   p = float(p)

   return int(0 - (n*math.log(p))/(math.log(2)**2))

def opt_ht(m, n):
   m = float(m)
   n = float(n)

   k = (m / n) * math.log(2)

   return int(max(1, round(k)))

def generate_read(n):
   read_list = []

   for i in range(n):
      read_list.append(random.choice(bases))

   return ''.join(read_list)

def get_neighbors(kmer_hash, K):
   neighbors = []
   kmer = khmer.reverse_hash(kmer_hash, K)

   begin = kmer[0:len(kmer)-1]
   end = kmer[1:len(kmer)]

   for base in bases:
      neighbors.append(khmer.forward_hash(base + begin, K))
      neighbors.append(khmer.forward_hash(end + base, K))

   return set(neighbors)

def explore(ht, start_kmer, K):
   discovered = set()
   explored = set()

   start_kmer_hash = khmer.forward_hash(start_kmer, K)
   
   if ht.get(kmer):
      discovered.add(start_kmer_hash)
   else:
      return 0

   while(len(discovered) > 0 and (len(explored) < 2000000)):
      kmer_hash = discovered.pop()
      kmer_neighbors = get_neighbors(kmer_hash, K)

      explored.add(kmer_hash)

      for neigh_hash in kmer_neighbors:
         if ht.get(khmer.reverse_hash(neigh_hash, K)) and neigh_hash not in explored and neigh_hash not in discovered:
            discovered.add(neigh_hash)

   return len(explored)

print "\"FPR\",\"LOWER\",\"AVG\",\"UPPER\""

K = 31
ps = [x/100.0 for x in range(1, 18)]
ps.append(0.1725)
for p in ps:
   sizes = []
   for i in range(10000):
      n = 1000
      contig_size = K
      m = calc_m(n, p)
      k = opt_ht(m, n)
      HT_SIZE = calc_ht_size(m, k)

      ht = khmer.new_hashbits(K, HT_SIZE, k)

      for j in range(n):
         kmer = generate_read(contig_size)
         ht.consume(kmer)

      kmer = generate_read(K)
      sizes.append(explore(ht, kmer, K))

   #avg = numpy.mean(sizes)
   #se = numpy.std(sizes) / numpy.sqrt(len(sizes))
   #lim = se * 1.96
   #print str(p) + "," + str(avg-lim) + "," + str(avg) + "," + str(avg+lim)
   low, med, upp = estimate_mean(sizes)
   print str(p) + "," + str(low) + "," + str(med) + "," + str(upp)
