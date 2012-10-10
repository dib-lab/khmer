import khmer
import screed
import sys
import math
import gc

K = 17

bases = ['A', 'C', 'G', 'T']

def callback(a, b, c):
   pass

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

def get_neighbors(kmer):
   neighbors = []

   begin = kmer[0:len(kmer)-1]
   end = kmer[1:len(kmer)]

   for base in bases:
      neighbors.append(base + begin)
      neighbors.append(end + base)

   return set(neighbors)

def find_neighbors(kmer, ht):
   tmp_neighs = get_neighbors(kmer)
   neighs = []

   for neigh in tmp_neighs:
      if ht.get(neigh):
         neighs.append(neigh)

   return neighs

def add_deg(degs, deg):
   if deg in degs.keys():
      degs[deg] += 1
   else:
      degs[deg] = 1

   return degs

def get_all_kmers(ht, start_kmer, K, ht2, degs):
   q = list()

   start_kmer_hash = khmer.forward_hash(start_kmer, K)

   if not ht2.get(start_kmer_hash):
      ht2.count(start_kmer)
   else:
      return ht2, degs

   neighs = find_neighbors(start_kmer, ht)

   degs = add_deg(degs, len(neighs))

   for neigh in neighs:
      neigh_hash = khmer.forward_hash(neigh, K)
      if not ht2.get(neigh):
         q.append(neigh_hash)
         ht2.count(neigh)

   counter = 0

   while len(q) != 0:
      counter += 1

      kmer_hash = q.pop()
      kmer = khmer.reverse_hash(kmer_hash, K)
      neighs = find_neighbors(kmer, ht)

      degs = add_deg(degs, len(neighs))

      for neigh in neighs:
         neigh_hash = khmer.forward_hash(neigh, K)
         if not ht2.get(neigh):
            q.append(neigh_hash)
            ht2.count(neigh)
         
   #return ht2, degs

def consume_and_save_ht(ht_size, n_ht, filename, outfile):
   ht = khmer.new_hashbits(K, ht_size, n_ht)
   ht.consume_fasta(filename, 0, 0, None, False, callback)
   ht.save(outfile + ".ht")
   return ht

def load_ht(filename):
   ht = khmer.new_hashbits(K, 1, 1)
   ht.load(filename)
   return ht

def diff(ht, filename):
   genome = khmer.new_hashbits(K, 4**K, 1)
   
   found = 0
   not_found = 0

   for n, record in enumerate(screed.fasta.fasta_iter(open(filename))):
      read = record['sequence']
      name = record['name']

      if 'N' in read:
         continue

      if len(read) < K:
         continue

      seq_len = len(read)
      for n in range(0,seq_len+1-K):
         kmer = read[n:n+K]

         if not genome.get(kmer):
            genome.consume(kmer)
         
            if ht.get(kmer):
               found += 1
            else:
               not_found += 1

   return found, not_found


def count(filename):
   ht = khmer.new_hashbits(K, 4**K, 1)
   ht2 = khmer.new_hashbits(K, 4**K, 1)
   ht3 = khmer.new_hashbits(K, 4**K, 1)

   unique_count = 0
   duplicate_count = 0
   neighbor_count = []
   
   for n, record in enumerate(screed.fasta.fasta_iter(open(filename))):
      read = record['sequence']
      name = record['name']

      if 'N' in read:
         continue

      if len(read) < K:
         continue

      seq_len = len(read)
      for n in range(0,seq_len+1-K):
         kmer = read[n:n+K]

         if ht.get(kmer) == 0:
            unique_count += 1
            ht.count(kmer)
         else:
            if ht2.get(kmer) == 0:
               duplicate_count += 1
               ht2.count(kmer)

   intercount = 0

   for n, record in enumerate(screed.fasta.fasta_iter(open(filename))):
      read = record['sequence']
      name = record['name']

      if len(read) < K:
         continue

      if 'N' in read:
         continue

      seq_len = len(read)
      for n in range(0,seq_len+1-K):
         kmer = read[n:n+K]

         if ht3.get(kmer) == 0:
            neighs = get_neighbors(kmer)
            ht3.count(kmer)

            neigh_count = 0

            for neigh in neighs:
               if ht.get(neigh) == 1:
                  neigh_count += 1
      
            if neigh_count > 2:
               intercount += 1
            neighbor_count.append(neigh_count)

   avg_deg = sum(neighbor_count) / float(len(neighbor_count))

   del ht
   del ht2
   del ht3

   return unique_count, duplicate_count, intercount, avg_deg

def deg(filename, ht):
   kmers = khmer.new_hashbits(K, 4**K, 1)

   degs = {}

   for n, record in enumerate(screed.fasta.fasta_iter(open(filename))):
      read = record['sequence']
      name = record['name']

      if len(read) < K:
         continue

      if 'N' in read:
         continue

      get_all_kmers(ht, read[0:K], K, kmers, degs)
      n_occ = kmers.n_occupied()
   
   del kmers

   return n_occ, degs

###

genome_file = sys.argv[1]
reads_file = sys.argv[2]

#ps = [0.01]
ps = [0.01, 0.05, 0.15]

#kmers = 4530123

kmers, dup_kmers, deg_kmers, avg_deg = count(genome_file)

print '%s\t%d\t%d\t%d\t%f\t%d\t%d' % ("ecoli 0%", kmers, 0, 0, 100, deg_kmers, (4**K)/8)

for p in ps:
   n = kmers
   m = calc_m(n, p)
   n_ht = opt_ht(m, n)
   ht_size = calc_ht_size(m, n_ht)
   outfile = str(p)

   ht = consume_and_save_ht(ht_size, n_ht, genome_file, outfile)
   found, not_found = diff(ht, genome_file)
   #print found, not_found, n, m, n_ht, ht_size
   num_found, deg_info = deg(genome_file, ht)

   deg_count = 0
   for i in range(3, 9):
      if i in deg_info.keys():
         deg_count += deg_info[i]

   num_addl = num_found - kmers + not_found
   perc_real = (kmers - not_found) / float(num_found)

   print '%s %f\t%d\t%d\t%d\t%f\t%d\t%d' % ("ecoli", p, num_found, num_addl, not_found, perc_real, deg_count, (ht_size*n_ht)/8)

   del ht
   gc.collect()


read_kmers, dup_read_kmers, deg_read_kmers, avg_deg_read = count(reads_file)
outfile = "readsexact.ht"
ht = consume_and_save_ht(4**K, 1, reads_file, outfile)
found, not_found = diff(ht, genome_file)
num_addl = read_kmers - kmers + not_found
perc_real = (kmers - not_found) / float(read_kmers)

print '%s\t%d\t%d\t%d\t%f\t%d\t%d' % ("reads 0%", read_kmers, num_addl, not_found, perc_real, deg_read_kmers, (4**K)/8)

del ht
gc.collect()

for p in ps:
   n = read_kmers
   m = calc_m(n, p)
   n_ht = opt_ht(m, n)
   ht_size = calc_ht_size(m, n_ht)
   outfile = "reads" + str(p)

   ht = consume_and_save_ht(ht_size, n_ht, reads_file, outfile)
   found, not_found = diff(ht, genome_file)
   num_found, deg_info = deg(reads_file, ht)

   deg_count = 0
   for i in range(3, 9):
      if i in deg_info.keys():
         deg_count += deg_info[i]

   num_addl = num_found - kmers + not_found
   perc_real = (kmers - not_found) / float(num_found)
   print '%s %f\t%d\t%d\t%d\t%f\t%d\t%d' % ("reads", p, num_found, num_addl, not_found, 
perc_real, deg_count, (ht_size*n_ht)/8)

   del ht
   gc.collect()
