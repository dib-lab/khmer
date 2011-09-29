__version__ = "0.2"

import _khmer
from _khmer import new_ktable
#from _khmer import new_hashtable
from _khmer import _new_counting_hash as new_hashtable
from _khmer import _new_counting_hash
from _khmer import _new_hashbits
from _khmer import new_readmask
from _khmer import new_minmax
#from _khmer import consume_genome
from _khmer import forward_hash, forward_hash_no_rc, reverse_hash
from _khmer import set_reporting_callback

from filter_utils import filter_fasta_file_any, filter_fasta_file_all, filter_fasta_file_limit_n

###

def new_hashbits(k, starting_size, n_tables=2):
    primes = get_n_primes_above_x(n_tables, starting_size)
    
    return _new_hashbits(k, primes)

def new_counting_hash(k, starting_size, n_tables=2):
    primes = get_n_primes_above_x(n_tables, starting_size)
    
    return _new_counting_hash(k, primes)

def load_hashbits(filename):
    ht = _new_hashbits(1, [1])
    ht.load(filename)

    return ht

def load_counting_hash(filename):
    ht = _new_counting_hash(1, [1])
    ht.load(filename)
    
    return ht

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
