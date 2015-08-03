#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2010-2015. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
"""This is khmer; please see http://khmer.readthedocs.org/."""

from __future__ import print_function
from math import log
import json

from khmer._khmer import CountingHash as _CountingHash
from khmer._khmer import LabelHash as _LabelHash
from khmer._khmer import Hashbits as _Hashbits
from khmer._khmer import HLLCounter as _HLLCounter
from khmer._khmer import ReadAligner as _ReadAligner

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
    hashtable = _CountingHash(1, [1])
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
    signature = None
    version = None
    ht_type = None

    uint_size = len(pack('I', 0))
    uchar_size = len(pack('B', 0))
    ulonglong_size = len(pack('Q', 0))

    try:
        with open(filename, 'rb') as hashbits:
            signature, = unpack('4s', hashbits.read(4))
            version, = unpack('B', hashbits.read(1))
            ht_type, = unpack('B', hashbits.read(1))
            ksize, = unpack('I', hashbits.read(uint_size))
            n_tables, = unpack('B', hashbits.read(uchar_size))
            table_size, = unpack('Q', hashbits.read(ulonglong_size))
        if signature != b"OXLI":
            raise ValueError("Node graph '{}' is missing file type "
                             "signature".format(filename) + str(signature))
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
    signature = None
    version = None
    ht_type = None
    use_bigcount = None

    uint_size = len(pack('I', 0))
    ulonglong_size = len(pack('Q', 0))

    try:
        with open(filename, 'rb') as countinghash:
            signature, = unpack('4s', countinghash.read(4))
            version, = unpack('B', countinghash.read(1))
            ht_type, = unpack('B', countinghash.read(1))
            use_bigcount, = unpack('B', countinghash.read(1))
            ksize, = unpack('I', countinghash.read(uint_size))
            n_tables, = unpack('B', countinghash.read(1))
            table_size, = unpack('Q', countinghash.read(ulonglong_size))
        if signature != b'OXLI':
            raise ValueError("Counting table '{}' is missing file type "
                             "signature. ".format(filename) + str(signature))
    except:
        raise ValueError("Counting table '{}' is corrupt ".format(filename))

    return ksize, round(table_size, -2), n_tables, use_bigcount, version, \
        ht_type


def calc_expected_collisions(hashtable, force=False, max_false_pos=.2):
    """Do a quick & dirty expected collision rate calculation on a hashtable.

    Also check to see that collision rate is within threshold.

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
        print("** ERROR: the graph structure is too small for ",
              file=sys.stderr)
        print("** this data set.  Increase data structure size",
              file=sys.stderr)
        print("** with --max_memory_usage/-M.", file=sys.stderr)
        print("**", file=sys.stderr)
        print("** Do not use these results!!", file=sys.stderr)
        print("**", file=sys.stderr)
        print("** (estimated false positive rate of %.3f;" % fp_all,
              file=sys.stderr, end=' ')
        print("max recommended %.3f)" % max_false_pos, file=sys.stderr)
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
    if target == 1 and number == 1:
        return [1]

    primes = []
    i = target - 1
    if i % 2 == 0:
        i -= 1
    while len(primes) != number and i > 0:
        if is_prime(i):
            primes.append(i)
        i -= 2

    if len(primes) != number:
        raise RuntimeError("unable to find %d prime numbers < %d" % (number,
                                                                     target))

    return primes


# Expose the cpython objects with __new__ implementations.
# These constructors add the functionality provided by the existing
# factory methods to the constructors defined over in cpython land.
# Additional functionality can be added to these classes as appropriate.


class CountingHash(_CountingHash):

    def __new__(cls, k, starting_size, n_tables):
        primes = get_n_primes_near_x(n_tables, starting_size)
        c = _CountingHash.__new__(cls, k, primes)
        c.primes = primes
        return c


class LabelHash(_LabelHash):

    def __new__(cls, k, starting_size, n_tables):
        hb = Hashbits(k, starting_size, n_tables)
        c = _LabelHash.__new__(cls, hb)
        c.graph = hb
        return c


class CountingLabelHash(_LabelHash):

    def __new__(cls, k, starting_size, n_tables):
        primes = get_n_primes_near_x(n_tables, starting_size)
        hb = _CountingHash(k, primes)
        c = _LabelHash.__new__(cls, hb)
        c.graph = hb
        return c


class Hashbits(_Hashbits):

    def __new__(cls, k, starting_size, n_tables):
        primes = get_n_primes_near_x(n_tables, starting_size)
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
        """Return the current estimation."""
        return _HLLCounter.estimate_cardinality(self)

    def estimate_cardinality(self):
        """Return the current estimation."""
        return _HLLCounter.estimate_cardinality(self)

    def add(self, kmer):
        """Add a k-mer to the counter."""
        return _HLLCounter.add(self, kmer)

    def consume_string(self, string):
        """Break a sequence into k-mers and add each k-mer to the counter."""
        return _HLLCounter.consume_string(self, string)

    def consume_fasta(self, fasta):
        """
        Reads sequences from file, breaks into kmers and adds kmers to counter.
        """
        return _HLLCounter.consume_fasta(self, fasta)

    def merge(self, other_counter):
        """Merge the other counter into this one."""
        return _HLLCounter.merge(self, other_counter)


class ReadAligner(_ReadAligner):

    """Sequence to graph aligner.

    ReadAligner uses a CountingHash (the counts of k-mers in the target DNA
    sequences) as an implicit De Bruijn graph. Input DNA sequences are aligned
    to this graph via a paired Hidden Markov Model.

    The HMM is configured upon class instantiation; default paramaters for the
    HMM are provided in 'defaultTransitionProbablitites' and
    'defaultScoringMatrix'.

    The main method is 'align'.
    """

    defaultTransitionProbabilities = (  # _M, _Ir, _Ig, _Mu, _Iru, _Igu
        (log(0.9848843, 2), log(0.0000735, 2), log(0.0000334, 2),
         log(0.0150068, 2), log(0.0000017, 2), log(0.0000003, 2)),  # M_
        (log(0.5196194, 2), log(0.4647955, 2), log(0.0059060, 2),
         log(0.0096792, 2)),  # Ir_
        (log(0.7611255, 2), log(0.2294619, 2), log(0.0072673, 2),
         log(0.0021453, 2)),  # Ig_
        (log(0.0799009, 2), log(0.0000262, 2), log(0.0001836, 2),
         log(0.9161349, 2), log(0.0033370, 2), log(0.0004173, 2)),  # Mu_
        (log(0.1434529, 2), log(0.0036995, 2), log(0.2642928, 2),
         log(0.5885548, 2)),  # Iru_
        (log(0.1384551, 2), log(0.0431328, 2), log(0.6362921, 2),
         log(0.1821200, 2))  # Igu_
    )

    defaultScoringMatrix = [
        log(0.955, 2), log(0.04, 2), log(0.004, 2), log(0.001, 2)]

    def __new__(cls, counting_table, trusted_cov_cutoff, bits_theta,
                **kwargs):

        if 'filename' in kwargs:
            with open(kwargs.pop('filename')) as paramfile:
                params = json.load(paramfile)
            scoring_matrix = params['scoring_matrix']
            transition_probabilities = params['transition_probabilities']
        else:
            if 'scoring_matrix' in kwargs:
                scoring_matrix = kwargs.pop('scoring_matrix')
            else:
                scoring_matrix = ReadAligner.defaultScoringMatrix
            if 'transition_probabilities' in kwargs:
                transition_probabilities = kwargs.pop(
                    'transition_probabilities')
            else:
                transition_probabilities = \
                    ReadAligner.defaultTransitionProbabilities
        r = _ReadAligner.__new__(cls, counting_table, trusted_cov_cutoff,
                                 bits_theta, scoring_matrix,
                                 transition_probabilities)
        r.graph = counting_table
        return r

    def __init__(self, *args, **kwargs):
        """
        ReadAligner initialization.

        HMM state notation abbreviations:
        M_t - trusted match; M_u - untrusted match
        Ir_t - trusted read insert; Ir_u - untrusted read insert
        Ig_t - trusted graph insert; Ig_u - untrusted graph insert

        Keyword arguments:
        filename - a path to a JSON encoded file providing the scoring matrix
            for the HMM in an entry named 'scoring_matrix' and the transition
            probababilties for the HMM in an entry named
            'transition_probabilities'. If provided the remaining keyword
            arguments are ignored. (default: None)
        scoring_matrix - a list of floats: trusted match, trusted mismatch,
            unstrusted match, untrusted mismatch. (default:
                ReadAligner.defaultScoringMatrix)
        transition_probabilities - A sparse matrix as a tuple of six tuples.
            The inner tuples contain 6, 4, 4, 6, 4, and 4 floats respectively.
            Transition are notated as 'StartState-NextState':
            (
              ( M_t-M_t,  M_t-Ir_t,  M_t-Ig_t,  M_t-M_u,  M_t-Ir_u,  M_t-Ig_u),
              (Ir_t-M_t, Ir_t-Ir_t,            Ir_t-M_u, Ir_t-Ir_u           ),
              (Ig_t-M_t,          , Ig_t-Ig_t, Ig_t-M_u,            Ig_t-Ig_u),
              ( M_u-M_t,  M_u-Ir_t,  M_u-Ig_t,  M_u-M_u,  M_u-Ir_u,  M_u-Ig_u),
              (Ir_u-M_t, Ir_u-Ir_t,            Ir_u-M_u, Ir_u-Ir_u           ),
              (Ig_u-M_t,          , Ig_u-Ig_t, Ig_u-M_u,            Ig_u-Ig_u)
            )
            (default: ReadAligner.defaultTransitionProbabilities)


        Note: the underlying CPython implementation creates the ReadAligner
        during the __new__ process and so the class initialization actually
        occurs there. Instatiation is documented here in __init__ as this is
        the traditional way.
        """
        _ReadAligner.__init__(self)

    def align_forward(self, read):
        return _ReadAligner.align_forward(self, read)

    def align(self, read):
        return _ReadAligner.align(self, read)

    def get_scoring_matrix(self):
        """
        Get the scoring matrix in use.

        Returns a tuple of floats:
        (trusted_match, trusted_mismatch, untrusted_match,untrusted_mismatch)
        """
        return _ReadAligner.get_scoring_matrix(self)

    def get_transition_probabilities(self):
        """
        Get the transition probabilties in use.

        HMM state notation abbreviations:
        M_t - trusted match;            M_u - untrusted match
        Ir_t - trusted read insert;     Ir_u - untrusted read insert
        Ig_t - trusted graph insert;    Ig_u - untrusted graph insert

        Returns a sparse matrix as a tuple of six tuples.
        The inner tuples contain 6, 4, 4, 6, 4, and 4 floats respectively.
        Transition are notated as 'StartState-NextState':

        {
            ( M_t-M_t,  M_t-Ir_t,  M_t-Ig_t,  M_t-M_u,  M_t-Ir_u,  M_t-Ig_u),
            (Ir_t-M_t, Ir_t-Ir_t,            Ir_t-M_u, Ir_t-Ir_u           ),
            (Ig_t-M_t,          , Ig_t-Ig_t, Ig_t-M_u,            Ig_t-Ig_u),
            ( M_u-M_t,  M_u-Ir_t,  M_u-Ig_t,  M_u-M_u,  M_u-Ir_u,  M_u-Ig_u),
            (Ir_u-M_t, Ir_u-Ir_t,            Ir_u-M_u, Ir_u-Ir_u           ),
            (Ig_u-M_t,          , Ig_u-Ig_t, Ig_u-M_u,            Ig_u-Ig_u)
        }
        """
        return _ReadAligner.get_transition_probabilities(self)
