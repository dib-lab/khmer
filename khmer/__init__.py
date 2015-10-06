# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2010-2015, Michigan State University.
# Copyright (C) 2015, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org
"""This is khmer; please see http://khmer.readthedocs.org/."""


from __future__ import print_function
from math import log
import json

from khmer._khmer import Countgraph as _Countgraph
from khmer._khmer import GraphLabels as _GraphLabels
from khmer._khmer import Nodegraph as _Nodegraph
from khmer._khmer import HLLCounter as _HLLCounter
from khmer._khmer import ReadAligner as _ReadAligner

from khmer._khmer import forward_hash
# tests/test_{functions,countgraph,counting_single}.py

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


def load_nodegraph(filename):
    """Load a nodegraph object from the given filename and return it.

    Keyword argument:
    filename -- the name of the nodegraph file
    """
    nodegraph = _Nodegraph(1, [1])
    nodegraph.load(filename)

    return nodegraph


def load_countgraph(filename):
    """Load a countgraph object from the given filename and return it.

    Keyword argument:
    filename -- the name of the countgraph file
    """
    countgraph = _Countgraph(1, [1])
    countgraph.load(filename)

    return countgraph


def extract_nodegraph_info(filename):
    """Open the given nodegraph file and return a tuple of information.

    Returns: the k-mer size, the table size, the number of tables, the version
    of the table format, and the type of table flag.

    Keyword argument:
    filename -- the name of the nodegraph file to inspect
    """
    ksize = None
    n_tables = None
    table_size = None
    signature = None
    version = None
    ht_type = None
    occupied = None

    uint_size = len(pack('I', 0))
    uchar_size = len(pack('B', 0))
    ulonglong_size = len(pack('Q', 0))

    try:
        with open(filename, 'rb') as nodegraph:
            signature, = unpack('4s', nodegraph.read(4))
            version, = unpack('B', nodegraph.read(1))
            ht_type, = unpack('B', nodegraph.read(1))
            ksize, = unpack('I', nodegraph.read(uint_size))
            n_tables, = unpack('B', nodegraph.read(uchar_size))
            occupied, = unpack('Q', nodegraph.read(ulonglong_size))
            table_size, = unpack('Q', nodegraph.read(ulonglong_size))
        if signature != b"OXLI":
            raise ValueError("Node graph '{}' is missing file type "
                             "signature".format(filename) + str(signature))
    except:
        raise ValueError("Node graph '{}' is corrupt ".format(filename))

    return ksize, round(table_size, -2), n_tables, version, ht_type, occupied


def extract_countgraph_info(filename):
    """Open the given countgraph file and return a tuple of information.

    Return: the k-mer size, the table size, the number of tables, the bigcount
    flag, the version of the table format, and the type of table flag.

    Keyword argument:
    filename -- the name of the countgraph file to inspect
    """
    ksize = None
    n_tables = None
    table_size = None
    signature = None
    version = None
    ht_type = None
    use_bigcount = None
    occupied = None

    uint_size = len(pack('I', 0))
    ulonglong_size = len(pack('Q', 0))

    try:
        with open(filename, 'rb') as countgraph:
            signature, = unpack('4s', countgraph.read(4))
            version, = unpack('B', countgraph.read(1))
            ht_type, = unpack('B', countgraph.read(1))
            use_bigcount, = unpack('B', countgraph.read(1))
            ksize, = unpack('I', countgraph.read(uint_size))
            n_tables, = unpack('B', countgraph.read(1))
            occupied, = unpack('Q', countgraph.read(ulonglong_size))
            table_size, = unpack('Q', countgraph.read(ulonglong_size))
        if signature != b'OXLI':
            raise ValueError("Count graph file '{}' is missing file type "
                             "signature. ".format(filename) + str(signature))
    except:
        raise ValueError("Count graph file '{}' is corrupt ".format(filename))

    return ksize, round(table_size, -2), n_tables, use_bigcount, version, \
        ht_type, occupied


def calc_expected_collisions(graph, force=False, max_false_pos=.2):
    """Do a quick & dirty expected collision rate calculation on a graph

    Also check to see that collision rate is within threshold.

    Keyword argument:
    graph: the countgraph or nodegraph object to inspect
    """
    sizes = graph.hashsizes()
    n_ht = float(len(sizes))
    occupancy = float(graph.n_occupied())
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


class Countgraph(_Countgraph):

    def __new__(cls, k, starting_size, n_tables):
        primes = get_n_primes_near_x(n_tables, starting_size)
        c = _Countgraph.__new__(cls, k, primes)
        c.primes = primes
        return c


class GraphLabels(_GraphLabels):

    def __new__(cls, k, starting_size, n_tables):
        hb = Nodegraph(k, starting_size, n_tables)
        c = _GraphLabels.__new__(cls, hb)
        c.graph = hb
        return c


class CountingGraphLabels(_GraphLabels):

    def __new__(cls, k, starting_size, n_tables):
        primes = get_n_primes_near_x(n_tables, starting_size)
        hb = _Countgraph(k, primes)
        c = _GraphLabels.__new__(cls, hb)
        c.graph = hb
        return c


class Nodegraph(_Nodegraph):

    def __new__(cls, k, starting_size, n_tables):
        primes = get_n_primes_near_x(n_tables, starting_size)
        c = _Nodegraph.__new__(cls, k, primes)
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


class ReadAligner(_ReadAligner):

    """Sequence to graph aligner.

    ReadAligner uses a Countgraph (the counts of k-mers in the target DNA
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

    def __new__(cls, count_graph, trusted_cov_cutoff, bits_theta,
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
        r = _ReadAligner.__new__(cls, count_graph, trusted_cov_cutoff,
                                 bits_theta, scoring_matrix,
                                 transition_probabilities)
        r.graph = count_graph
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
