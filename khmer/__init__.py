#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#

from _khmer import new_ktable
from _khmer import _new_counting_hash
from _khmer import _new_hashbits
from _khmer import set_reporting_callback
from _khmer import new_readaligner
from _khmer import forward_hash
from _khmer import new_hashtable
from _khmer import forward_hash_no_rc
from _khmer import reverse_hash
from _khmer import get_config
from _khmer import ReadParser
from _khmer import _LabelHash
from _khmer import _Hashbits

from struct import pack, unpack

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

#


def new_hashbits(k, starting_size, n_tables=2):
    primes = get_n_primes_above_x(n_tables, starting_size)

    return _new_hashbits(k, primes)


def new_counting_hash(k, starting_size, n_tables=2, n_threads=1):
    primes = get_n_primes_above_x(n_tables, starting_size)

    return _new_counting_hash(k, primes, n_threads)


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


def extract_hashbits_info(filename):
    ksize = None
    n_tables = None
    table_size = None
    version = None
    ht_type = None

    uint_size = len(pack('I', 0))
    uchar_size = len(pack('B', 0))
    ulonglong_size = len(pack('Q', 0))

    with open(filename, 'rb') as f:
        version, = unpack('B', f.read(1))
        ht_type, = unpack('B', f.read(1))
        ksize, = unpack('I', f.read(uint_size))
        n_tables, = unpack('B', f.read(uchar_size))
        table_size, = unpack('Q', f.read(ulonglong_size))

    return ksize, round(table_size, -2), n_tables, version, ht_type


def extract_countinghash_info(filename):
    ksize = None
    n_tables = None
    table_size = None
    version = None
    ht_type = None
    use_bigcount = None

    uint_size = len(pack('I', 0))
    uchar_size = len(pack('B', 0))
    ulonglong_size = len(pack('Q', 0))

    with open(filename, 'rb') as f:
        version, = unpack('B', f.read(1))
        ht_type, = unpack('B', f.read(1))
        use_bigcount, = unpack('B', f.read(1))
        ksize, = unpack('I', f.read(uint_size))
        n_tables, = unpack('B', f.read(1))
        table_size, = unpack('Q', f.read(ulonglong_size))

    return ksize, round(table_size, -2), n_tables, \
        use_bigcount, version, ht_type


def calc_expected_collisions(ht):
    """
    A quick & dirty expected collision rate calculation
    """
    sizes = ht.hashsizes()
    n_ht = float(len(sizes))
    occupancy = float(ht.n_occupied())
    min_size = min(sizes)

    fp_one = occupancy / min_size
    fp_all = fp_one ** n_ht

    return fp_all

#


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
    for x in range(3, int(n ** 0.5) + 1, 2):
        if n % x == 0:
            return False
    return True


def get_n_primes_near_x(n, x):
    '''
    steps backward until n primes (other than 2) have been
    found that are smaller than x.
    '''
    primes = []
    i = x - 1
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
    i = x + 1
    if i % 2 == 0:
        i += 1
    while len(primes) != n and i > 0:
        if is_prime(i):
            primes.append(i)
        i += 2
    return primes

'''
Expose the cpython objects with __new__ implementations.
These constructors add the functionality provided by the existing
factory methods to the constructors defined over in cpython land.
Additional functionality can be added to these classes as appropriate.
'''


class LabelHash(_LabelHash):

    def __new__(cls, k, starting_size, n_tables):
        primes = get_n_primes_above_x(n_tables, starting_size)
        c = _LabelHash.__new__(cls, k, primes)
        c.primes = primes
        return c


class Hashbits(_Hashbits):

    def __new__(cls, k, starting_size, n_tables):
        primes = get_n_primes_above_x(n_tables, starting_size)
        c = _Hashbits.__new__(cls, k, primes)
        c.primes = primes
        return c
