__version__ = "0.2"

import _khmer
from _khmer import new_ktable
from _khmer import new_hashtable
from _khmer import consume_genome

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

# from http://www.rsok.com/~jrm/printprimes.html
PRIMES_1m = [1000003, 1009837]
PRIMES_100m = [100009979, 100000007]

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
