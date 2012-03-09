import khmer
import sys
import random
import math

def callback(a, b, c):
   pass

bases = ['A', 'C', 'G', 'T']

#p = float(sys.argv[1])

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
   edges = set()
   discovered = set()
   explored = set()
   hash_ids = {}

   start_kmer_hash = khmer.forward_hash(start_kmer, K)

   if ht.get(khmer.reverse_hash(start_kmer_hash, K)):
      discovered.add(start_kmer_hash)
      hash_ids[start_kmer_hash] = len(hash_ids.keys()) + 1
   else:
      return hash_ids, edges

   while(len(discovered) > 0 and (len(explored) < 2000000)):
      kmer_hash = discovered.pop()
      kmer_neighbors = get_neighbors(kmer_hash, K)

      explored.add(kmer_hash)

      for neigh_hash in kmer_neighbors:
         if ht.get(khmer.reverse_hash(neigh_hash, K)) and neigh_hash not in explored and neigh_hash not in discovered:
            discovered.add(neigh_hash)
            hash_ids[neigh_hash] = len(hash_ids.keys()) + 1
            edges.add(tuple(sorted([hash_ids[neigh_hash], hash_ids[kmer_hash]])))
         elif ht.get(khmer.reverse_hash(neigh_hash, K)) and (neigh_hash in explored or neigh_hash in discovered):
            edges.add(tuple(sorted([hash_ids[neigh_hash], hash_ids[kmer_hash]])))
         
   return hash_ids, edges 

def gen_graph(filename, edges, hash_ids, chr, K):
   fd = open(filename, "w")

   fd.write("graph x {\nsize=\"16, 16\";\n")
   fd.write("node [ color = red, fontcolor = black, style = filled ];\n")

   for i in range(len(chr) - K):
      kmer = chr[i:i + K]
      kmer_hash = khmer.forward_hash(kmer, K)
      hash_id = hash_ids[kmer_hash]

      fd.write("N" + str(hash_id) + " [color = black, fontcolor = white];\n")

   for edge in edges:
      fd.write("N" + str(edge[0]) + " -- " + "N" + str(edge[1]) + ";\n")

   fd.write("}")
   fd.close()

def main():
   K = 31
   n = 1000
   seq = gen_circ_chrom(n, K)

   for p in [0.01, 0.05, 0.1, 0.15]:
      m = calc_m(n, p)
      k = opt_ht(m, n)
      HT_SIZE = calc_ht_size(m, k)

      ht = khmer.new_hashbits(K, HT_SIZE, k)
      ht.consume(seq)
      hash_ids, edges = explore(ht, seq[0:K], K)

      gen_graph(str(p) + ".dot", edges, hash_ids, seq, K)

main()
