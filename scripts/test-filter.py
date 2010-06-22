import sys
import glob
import khmer

K = 17
THRESHOLD = 3
MINLENGTH = 50

def fib(n):
   if n == 0 or n == 1:
      return 1
   else:
      return fib(n-1) + fib(n-2)

hasher = khmer.new_hashtable(K, 16000000000)

fib(100)

#files = glob.glob(sys.argv[1])

#for filename in files:
#   hasher.consume_fasta(filename)

#for filename in files:
#   hasher.filter_fasta_file(filename, MINLENGTH, THRESHOLD);
