__version__ = "0.2"

import _khmer
from _khmer import new_ktable
from _khmer import new_hashtable
from _khmer import _new_hashbits
from _khmer import new_readmask
from _khmer import new_minmax
from _khmer import consume_genome
from _khmer import forward_hash, forward_hash_no_rc, reverse_hash
from _khmer import set_reporting_callback
from _khmer import do_intersection_partition

from filter_utils import filter_fasta_file_any, filter_fasta_file_all, filter_fasta_file_limit_n

###

def new_hashbits(k, starting_size, n_tables=8):
    primes = get_n_primes_above_x(n_tables, starting_size)
    print primes
#    primes = [ 22906493, 22906519, 22906561, 22906567,
#               22906619, 22906649, 22906657, 22906661 ]
    
    return _new_hashbits(k, primes)

def _default_reporting_callback(info, n_reads, other):
    print '...', info, n_reads, other

def reset_reporting_callback():
    set_reporting_callback(_default_reporting_callback)

reset_reporting_callback()

###

class KmerCount(object):
    def __init__(self, size, report_zero=False):
        self._kt = new_ktable(size)
        self.report_zero = report_zero

    def consume(self, seq):
        self._kt.consume(seq)

    def _get_pairs(self):
        kt = self._kt
        size = kt.n_entries()
        
        for i in range(0, size):
            count = kt.get(i)
            if count or self.report_zero:
                kmer = kt.reverse_hash(i)
                yield kmer, kt.get(i)

    pairs = property(_get_pairs)

    def __getitem__(self, k):
        return self._kt.get(k)

def is_prime(n):
   '''
   checks if a number is prime
   '''
   if n < 2:
      return False
   if n == 2:
      return True
   if n % 2 == 0:
      return False
   for x in range(3, int(n**0.5)+1, 2):
      if n % x == 0:
         return False
   return True

def get_n_primes_near_x(n, x):
   '''
   steps backward until n primes (other than 2) have been
   found that are smaller than x.
   '''
   primes = []
   i = x-1
   if i % 2 == 0:
      i -= 1
   while len(primes) != n and i > 0:
      if is_prime(i):
         primes.append(i)
      i -= 2
   return primes

def get_n_primes_above_x(n, x):
   '''
   steps forward until n primes (other than 2) have been
   found that are smaller than x.
   '''
   primes = []
   i = x+1
   if i % 2 == 0:
      i += 1
   while len(primes) != n and i > 0:
      if is_prime(i):
         primes.append(i)
      i += 2
   return primes

# from http://www.rsok.com/~jrm/printprimes.html
PRIMES_1m = [1000003, 1009837]
PRIMES_100m = [100009979, 100000007]
PRIMES_1b = [1000000007, 1000000919]
PRIMES_2b = [1999999973, 1999999943]
PRIMES_4b = [4000000007, 4000000009]
PRIMES_8b = [8000000011, 8000000051]

class HashtableIntersect(object):
    def __init__(self, k, size1, size2):
        self._kh1 = new_hashtable(k, size1)
        self._kh2 = new_hashtable(k, size2)

    def consume(self, seq):
        self._kh1.consume(seq)
        self._kh2.consume(seq)

    def get_min_count(self, seq):
        return min(self._kh1.get_min_count(seq),
                   self._kh2.get_min_count(seq))

    def get_max_count(self, seq):
        return min(self._kh1.get_max_count(seq),
                   self._kh2.get_max_count(seq))
