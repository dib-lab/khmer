import sys
import glob
import khmer

K = 17
THRESHOLD = 5
MINLENGTH = 50

infile = sys.argv[1]
#outfile = sys.argv[2]

primes = [8000000111, 8000000257, 8000000321, 8000000381, 8000000423]

hasher = khmer.new_hashtable(K, 8000000011)
hasher.consume_fasta(infile)
hasher.filter_fasta_file(infile, "8000000011.fa", MINLENGTH, THRESHOLD)

prevPrime = 8000000011

for prime in primes:
   hasher = khmer.new_hashtable(K, prime)

   hasher.consume_fasta(str(prevPrime) + ".fa")

   hasher.filter_fasta_file(str(prevPrime) + ".fa", str(prime) + ".fa", MINLENGTH, THRESHOLD);

   prevPrime = prime
