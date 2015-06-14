#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2010-2015. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
"""This is khmer; please see http://khmer.readthedocs.org/."""

from __future__ import print_function

from khmer._khmer import CountingHash
from khmer._khmer import LabelHash as _LabelHash
from khmer._khmer import Hashbits as _Hashbits
from khmer._khmer import HLLCounter as _HLLCounter
from khmer._khmer import ReadAligner

from khmer._khmer import forward_hash  # figuregen/*.py
# tests/test_{functions,counting_hash,labelhash,counting_single}.py

from khmer._khmer import new_hashtable
# sandbox/{occupy,ctb-iterative-bench{-2-old}}.py
# tests/{test_c_wrapper,test_counting_single}.py

from khmer._khmer import forward_hash_no_rc  # tests/test_functions.py

from khmer._khmer import reverse_hash  # tests/test_functions.py
# tests/counting_single.py

from khmer._khmer import hash_murmur3        # tests/test_functions.py
from khmer._khmer import hash_no_rc_murmur3  # tests/test_functions.py

from khmer._khmer import get_version_cpp as __version_cpp__
# tests/test_version.py

from khmer._khmer import ReadParser  # sandbox/to-casava-1.8-fastq.py
# tests/test_read_parsers.py,scripts/{filter-abund-single,load-graph}.py
# scripts/{abundance-dist-single,load-into-counting}.py

import sys

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

    return _Hashbits(k, primes)


def new_counting_hash(k, starting_size, n_tables=2):
    """Return a new countinghash object.

    Keyword arguments:
    k -- kmer size to use
    starting_size -- lower bound on hashsize to use
    n_tables -- number of hash tables to use (default = 2)
    n_threads  -- number of simultaneous threads to execute (default = 1)
    """
    primes = get_n_primes_above_x(n_tables, starting_size)

    return CountingHash(k, primes)


def load_hashbits(filename):
    """Load a hashbits object from the given filename and return it.

    Keyword argument:
    filename -- the name of the hashbits file
    """
    hashtable = _Hashbits(1, [1])
    hashtable.load(filename)

    return hashtable


def load_counting_hash(filename):
    """Load a counting_hash object from the given filename and return it.

    Keyword argument:
    filename -- the name of the counting_hash file
    """
    hashtable = CountingHash(1, [1])
    hashtable.load(filename)

    return hashtable


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

    try:
        with open(filename, 'rb') as hashbits:
            version, = unpack('B', hashbits.read(1))
            ht_type, = unpack('B', hashbits.read(1))
            ksize, = unpack('I', hashbits.read(uint_size))
            n_tables, = unpack('B', hashbits.read(uchar_size))
            table_size, = unpack('Q', hashbits.read(ulonglong_size))
    except:
        raise ValueError("Presence table '{}' is corrupt ".format(filename))

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

    try:
        with open(filename, 'rb') as countinghash:
            version, = unpack('B', countinghash.read(1))
            ht_type, = unpack('B', countinghash.read(1))
            use_bigcount, = unpack('B', countinghash.read(1))
            ksize, = unpack('I', countinghash.read(uint_size))
            n_tables, = unpack('B', countinghash.read(1))
            table_size, = unpack('Q', countinghash.read(ulonglong_size))
    except:
        raise ValueError("Counting table '{}' is corrupt ".format(filename))

    return ksize, round(table_size, -2), n_tables, use_bigcount, version, \
        ht_type


def calc_expected_collisions(hashtable, force=False, max_false_pos=.2):
    """Do a quick & dirty expected collision rate calculation on a hashtable.
    Check to see that collision rate is within threshold.

    Keyword argument:
    hashtable: the hashtable object to inspect
    """
    sizes = hashtable.hashsizes()
    n_ht = float(len(sizes))
    occupancy = float(hashtable.n_occupied())
    min_size = min(sizes)

    fp_one = occupancy / min_size
    fp_all = fp_one ** n_ht

    if fp_all > max_false_pos:
        print("**", file=sys.stderr)
        print(
            "** ERROR: the graph structure is too small for ", file=sys.stderr)
        print(
            "this data set.  Increase k-mer presence table ", file=sys.stderr)
        print("size/num of tables.", file=sys.stderr)
        print("** Do not use these results!!", file=sys.stderr)
        print("**", file=sys.stderr)
        if not force:
            sys.exit(1)

    return fp_all


def is_prime(number):
    """Check if a number is prime."""
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
    """Backward-find primes smaller than target.

    Step backwards until a number of primes (other than 2) have been
    found that are smaller than the target and return them.

    Keyword arguments:
    number -- the number of primes to find
    target -- the number to step backwards from
    """
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
    """Forward-find primes smaller than target.

    Step forwards until a number of primes (other than 2) have been
    found that are smaller than the target and return them.

    Keyword arguments:
    number -- the number of primes to find
    target -- the number to step forwards from
    """
    primes = []
    i = target + 1
    if i % 2 == 0:
        i += 1
    while len(primes) != number and i > 0:
        if is_prime(i):
            primes.append(i)
        i += 2
    return primes


# Expose the cpython objects with __new__ implementations.
# These constructors add the functionality provided by the existing
# factory methods to the constructors defined over in cpython land.
# Additional functionality can be added to these classes as appropriate.


class LabelHash(_LabelHash):

    def __new__(cls, k, starting_size, n_tables):
        hb = Hashbits(k, starting_size, n_tables)
        c = _LabelHash.__new__(cls, hb)
        c.graph = hb
        return c


class CountingLabelHash(_LabelHash):

    def __new__(cls, k, starting_size, n_tables):
        primes = get_n_primes_above_x(n_tables, starting_size)
        hb = CountingHash(k, primes)
        c = _LabelHash.__new__(cls, hb)
        c.graph = hb
        return c


class Hashbits(_Hashbits):

    def __new__(cls, k, starting_size, n_tables):
        primes = get_n_primes_above_x(n_tables, starting_size)
        c = _Hashbits.__new__(cls, k, primes)
        c.primes = primes
        return c


class HLLCounter(_HLLCounter):

    """HyperLogLog counter.

    A HyperLogLog counter is a probabilistic data structure specialized on
    cardinality estimation.
    There is a precision/memory consumption trade-off: error rate determines
    how much memory is consumed.

    # Creating a new HLLCounter:

    >>> khmer.HLLCounter(error_rate, ksize)

    where the default values are:
      - error_rate: 0.01
      - ksize: 20
    """

    def __len__(self):
        return self.estimate_cardinality()
