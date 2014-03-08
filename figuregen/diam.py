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
import copy
import numpy

from collections import deque

def callback(a, b, c):
   pass

bases = ['A', 'C', 'G', 'T']

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

def gen_circ_chrom(n, K):
   read = generate_read(n)
   chromosome = read + read[0:K]
   return chromosome

def get_neighbors(kmer, K):
   neighbors = []

   begin = kmer[0:len(kmer)-1]
   end = kmer[1:len(kmer)]

   for base in bases:
      neighbors.append(base + begin)
      neighbors.append(end + base)

   return set(neighbors)

def find_neighbors(kmer, K, ht):
   tmp_neighs = get_neighbors(kmer, K)
   neighs = []

   for neigh in tmp_neighs:
      if ht.get(neigh):
         neighs.append(neigh)

   return neighs

def get_all_kmers(ht, start_kmer, K):
   q = []
   verts = {}

   level = 0

   verts[start_kmer] = level

   neighs = find_neighbors(start_kmer, K, ht)

   for neigh in neighs:
         q.append(neigh)

   while len(q) != 0:
      new_q = []
      level += 1

      for kmer in q:
         neighs = find_neighbors(kmer, K, ht)
         for neigh in neighs:
            if neigh not in verts.keys():
               new_q.append(neigh)

         verts[kmer] = level

      q = new_q

   return verts.keys()

def get_real_kmers(seq, K):
   kmers = set()
   n = len(seq) - K + 1

   for i in range(n):
      kmer = seq[i:i+K]
      kmers.add(kmer)

   return kmers

def get_level(ht, start_kmer, real, K):
   q = []
   vert_set = set()

   real_kmers = copy.deepcopy(real)

   level = 0

   vert_set.add(start_kmer)
   if start_kmer in real_kmers:
      real_kmers.remove(start_kmer)
 
   neighs = find_neighbors(start_kmer, K, ht)

   for neigh in neighs:
      q.append(neigh)

   while len(q) != 0:
      new_q = []

      if len(real_kmers) == 0:
         return level

      level += 1

      for kmer in q:
         neighs = find_neighbors(kmer, K, ht)
         for neigh in neighs:
            if neigh not in vert_set:
               new_q.append(neigh)

         vert_set.add(kmer)
         if kmer in real_kmers:
            real_kmers.remove(kmer)

      q = new_q

   return level

def main():
   K = 8
   n = 50
   add_kmers = 50
   total_kmers = n + add_kmers

   print "\"FPR\",\"LOWER\",\"AVG\",\"UPPER\""

   for p in [x/200.0 + .01 for x in range(59)]:
      diam_lens = []
      for j in range(500):
         seq = gen_circ_chrom(n, K)
         m = calc_m(total_kmers, p)
         k = opt_ht(m, total_kmers)
         HT_SIZE = calc_ht_size(m, k)

         ht = khmer.new_hashbits(K, HT_SIZE, k)
         ht.consume(seq)

         for i in range(add_kmers):
            ht.consume(generate_read(K))

         real_kmers = get_real_kmers(seq, K)

         out_len = []
         # step one: find the "outbranch" lengths for each real k-mer
         for kmer in real_kmers:
            out_len.append(get_level(ht, kmer, real_kmers, K))

         # step two: find the shortest longest path using the info from step 1
         diam_lens.append(max(out_len))

      #avg = numpy.mean(diam_lens)
      #se = numpy.std(diam_lens) / numpy.sqrt(len(diam_lens))
      #lim = se * 1.96
      #print str(p) + "," + str(avg-lim) + "," + str(avg) + "," + str(avg+lim)
      low, med, upp = estimate_mean(diam_lens)
      print str(p) + "," + str(low) + "," + str(med) + "," + str(upp)

main()
