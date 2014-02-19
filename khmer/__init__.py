'''
This file is part of khmer, http://github.com/ged-lab/khmer/, and is
Copyright (C) Michigan State University, 2009-2014. It is licensed under
the three-clause BSD license; see doc/LICENSE.txt.
Contact: khmer-project@idyll.org
'''

from khmer._khmer import new_ktable
from khmer._khmer import _new_counting_hash
from khmer._khmer import _new_hashbits
from khmer._khmer import set_reporting_callback
from khmer._khmer import _LabelHash
from khmer._khmer import _Hashbits

from struct import pack, unpack

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


def new_hashbits(k, starting_size, n_tables=2):
    """Return a new hashbits object. Deprecated.

    This factory method is deprecated in favor of creating a Hashbits object
    directly via 'new Hashbits(...)'.

    Keyword argument:
    k -- kmer size to use
    starting_size -- lower bound on hashsize to use
    n_tables -- number of hash tables to use (default = 2)
    """
    primes = get_n_primes_above_x(n_tables, starting_size)

    return _new_hashbits(k, primes)


def new_counting_hash(k, starting_size, n_tables=2, n_threads=1):
    """Return a new countinghash object.

    Keyword arguments:
    k -- kmer size to use
    starting_size -- lower bound on hashsize to use
    n_tables -- number of hash tables to use (default = 2)
    n_threads  -- number of simultaneous threads to execute (default = 1)
    """
    primes = get_n_primes_above_x(n_tables, starting_size)

    return _new_counting_hash(k, primes, n_threads)


def load_hashbits(filename):
    """Load a hashbits object from the given filename and return it.

    Keyword argument:
    filename -- the name of the hashbits file
    """
    hashtable = _new_hashbits(1, [1])
    hashtable.load(filename)

    return hashtable


def load_counting_hash(filename):
    """Load a counting_hash object from the given filename and return it.

    Keyword argument:
    filename -- the name of the counting_hash file
    """
    hashtable = _new_counting_hash(1, [1])
    hashtable.load(filename)

    return hashtable


def _default_reporting_callback(info, n_reads, other):
    print '...', info, n_reads, other


def reset_reporting_callback():
    set_reporting_callback(_default_reporting_callback)

reset_reporting_callback()


def extract_hashbits_info(filename):
    """Open the given hashbits file and return a tuple of information.

    Returns: the k-mer size, the table size, the number of tables, the version
    of the table format, and the type of table flag.

    Keyword argument:
    filename -- the name of the hashbits file to inspect
    """
    ksize = None
    n_tables = None
    table_size = None
    version = None
    ht_type = None

    uint_size = len(pack('I', 0))
    uchar_size = len(pack('B', 0))
    ulonglong_size = len(pack('Q', 0))

    with open(filename, 'rb') as hashbits:
        version, = unpack('B', hashbits.read(1))
        ht_type, = unpack('B', hashbits.read(1))
        ksize, = unpack('I', hashbits.read(uint_size))
        n_tables, = unpack('B', hashbits.read(uchar_size))
        table_size, = unpack('Q', hashbits.read(ulonglong_size))

    return ksize, round(table_size, -2), n_tables, version, ht_type


def extract_countinghash_info(filename):
    """Open the given counting_hash file and return a tuple of information.

    Return: the k-mer size, the table size, the number of tables, the bigcount
    flag, the version of the table format, and the type of table flag.

    Keyword argument:
    filename -- the name of the counting_hash file to inspect
    """
    ksize = None
    n_tables = None
    table_size = None
    version = None
    ht_type = None
    use_bigcount = None

    uint_size = len(pack('I', 0))
    ulonglong_size = len(pack('Q', 0))

    with open(filename, 'rb') as countinghash:
        version, = unpack('B', countinghash.read(1))
        ht_type, = unpack('B', countinghash.read(1))
        use_bigcount, = unpack('B', countinghash.read(1))
        ksize, = unpack('I', countinghash.read(uint_size))
        n_tables, = unpack('B', countinghash.read(1))
        table_size, = unpack('Q', countinghash.read(ulonglong_size))

    return ksize, round(table_size, -2), n_tables, use_bigcount, version, \
        ht_type


def calc_expected_collisions(hashtable):
    """Do a quick & dirty expected collision rate calculation on a hashtable.

    Keyword argument:
    hashtable: the hashtable object to inspect
    """
    sizes = hashtable.hashsizes()
    n_ht = float(len(sizes))
    occupancy = float(hashtable.n_occupied())
    min_size = min(sizes)

    fp_one = occupancy / min_size
    fp_all = fp_one ** n_ht

    return fp_all


class KmerCount(object):

    def __init__(self, size, report_zero=False):
        self._kt = new_ktable(size)
        self.report_zero = report_zero

    def consume(self, seq):
        self._kt.consume(seq)

    def _get_pairs(self):
        ktable = self._kt
        size = ktable.n_entries()

        for i in range(0, size):
            count = ktable.get(i)
            if count or self.report_zero:
                kmer = ktable.reverse_hash(i)
                yield kmer, ktable.get(i)

    pairs = property(_get_pairs)

    def __getitem__(self, k):
        return self._kt.get(k)


def is_prime(number):
    '''Checks if a number is prime.'''
    if number < 2:
        return False
    if number == 2:
        return True
    if number % 2 == 0:
        return False
    for _ in range(3, int(number ** 0.5) + 1, 2):
        if number % _ == 0:
            return False
    return True


def get_n_primes_near_x(number, target):
    ''' Step backwards until a number of primes (other than 2) have been
    found that are smaller than the target and return them.

    Keyword arguments:
    number -- the number of primes to find
    target -- the number to step backwards from
    '''
    primes = []
    i = target - 1
    if i % 2 == 0:
        i -= 1
    while len(primes) != number and i > 0:
        if is_prime(i):
            primes.append(i)
        i -= 2
    return primes


def get_n_primes_above_x(number, target):
    '''Step forwards until a number of primes (other than 2) have been
    found that are smaller than the target and return them.

    Keyword arguments:
    number -- the number of primes to find
    target -- the number to step forwards from
    '''
    primes = []
    i = target + 1
    if i % 2 == 0:
        i += 1
    while len(primes) != number and i > 0:
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
