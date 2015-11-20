#! /usr/bin/env python

"""
Hypothesis tests for khmer.

Hypothesis is a library implementing a mix of property-based testing,
unit testing and fuzzing. Tests are similar to unit tests,
but specify strategies for generating input.
Hypothesis takes these strategies and tries to find inputs to falsify the tests.
"""

from __future__ import division, unicode_literals

from collections import Counter

from hypothesis import given, assume, strategies as st
from nose.plugins.attrib import attr

import khmer
from khmer import reverse_hash, forward_hash, forward_hash_no_rc


# TODO: we are only testing with fixed k, table size and number of tables for now.
KSIZE = 13
N_TABLES = 4
TABLE_SIZE = 4 ** 4

# strategy for creating kmers. Alphabet is derived from nucleotides.
st_kmer = st.text("ACGT", min_size=KSIZE, max_size=KSIZE)

# Reverse complement utilities.
TRANSLATE = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
def revcomp(kmer):
    return "".join(TRANSLATE[c] for c in kmer[::-1])


def lex_rc(kmer):
    """ Returns the smaller string between kmer and its reverse complement,
    in lexicographical order """
    # We abuse one property from reverse_hash here,
    # since its implementation always return the smaller between kmer and rc.
    return reverse_hash(forward_hash(kmer, len(kmer)), len(kmer))


@attr('hypothesis')
@given(st_kmer)
def test_forward_hash(kmer):
    ksize = len(kmer)

    rh = reverse_hash(forward_hash(kmer, ksize), ksize)

    # We need to check for either kmer or reverse complement,
    # due to how the hash functions are implemented.
    assert rh == kmer or rh == revcomp(kmer)


@attr('hypothesis')
@given(st_kmer)
def test_forward_hash_no_rc(kmer):
    ksize = len(kmer)

    rh = reverse_hash(forward_hash_no_rc(kmer, ksize), ksize)

    # forward_hash_no_rc only generates the hash for the kmer,
    # not the reverse complement. In this case we just need to check if
    # the kmer is the same.
    assert rh == kmer


@attr('hypothesis')
@given(st.lists(st_kmer, min_size=1))
def test_countgraph_undercounting(kmers):
    """Testing countgraph undercounting.

    A collections.Counter serves as an oracle for Count-Min sketches,
    since both implement a frequency counter interface."""

    oracle = Counter()
    countgraph = khmer.Countgraph(KSIZE, TABLE_SIZE, N_TABLES)

    for kmer in kmers:
        oracle.update([kmer])
        countgraph.count(kmer)

    for kmer in oracle:
        # Our implementation only counts to 255 by default,
        # so we need to check for at most 255 in the comparison.
        assert countgraph.get(kmer) >= min(oracle[kmer], 255)


@attr('hypothesis')
@given(st.lists(st_kmer, min_size=1))
def test_countgraph_undercounting_consume(kmers):
    """Testing countgraph undercounting, using consume instead of count.
        
        A collections.Counter serves as an oracle for Count-Min sketches,
        since both implement a frequency counter interface."""
    
    oracle = Counter()
    countgraph = khmer.Countgraph(KSIZE, TABLE_SIZE, N_TABLES)
    
    for kmer in kmers:
        oracle.update([kmer])
        countgraph.consume(kmer)
    
    for kmer in oracle:
        # Our implementation only counts to 255 by default,
        # so we need to check for at most 255 in the comparison.
        assert countgraph.get(kmer) >= min(oracle[kmer], 255)


@attr('hypothesis')
@given(st.lists(st_kmer, min_size=1))
def test_countgraph_undercounting_bigcounts(kmers):
    """Testing countgraph undercounting, using the bigcount feature.

    A collections.Counter serves as an oracle for Count-Min sketches,
    since both implement a frequency counter interface."""

    oracle = Counter()
    countgraph = khmer.Countgraph(KSIZE, TABLE_SIZE, N_TABLES)
    countgraph.set_use_bigcount(True)

    for kmer in kmers:
        oracle.update([kmer])
        countgraph.count(kmer)

    for kmer in oracle:
        # The bigcount implementation counts to 65535,
        # so we need to check for this at most in the comparison.
        assert countgraph.get(kmer) >= min(oracle[kmer], 65535)


@attr('hypothesis')
@given(st.lists(st_kmer, min_size=1))
def test_countgraph_undercounting_bigcounts_consume(kmers):
    """Testing countgraph undercounting, using the bigcount feature and using consume.

    A collections.Counter serves as an oracle for Count-Min sketches,
    since both implement a frequency counter interface."""

    oracle = Counter()
    countgraph = khmer.Countgraph(KSIZE, TABLE_SIZE, N_TABLES)
    countgraph.set_use_bigcount(True)

    for kmer in kmers:
        oracle.update([kmer])
        countgraph.consume(kmer)

    for kmer in oracle:
        # The bigcount implementation counts to 65535,
        # so we need to check for this at most in the comparison.
        assert countgraph.get(kmer) >= min(oracle[kmer], 65535)


@attr('hypothesis')
@given(st.sets(st_kmer, min_size=1))
def test_nodegraph_presence(kmers):
    """Testing nodegraph for presence checking.

    A set serves as an oracle for Bloom Filters,
    since both implement set membership operations."""

    oracle = set()
    nodegraph = khmer.Nodegraph(KSIZE, TABLE_SIZE, N_TABLES)

    for kmer in kmers:
        oracle.update([kmer])
        nodegraph.count(kmer)

    for kmer in oracle:
        # Nodegraph and Countgraph have similar API,
        # but Nodegraph.get returns only {0,1}
        assert nodegraph.get(kmer) == 1


@attr('hypothesis')
@attr('known_failing')
@given(st.lists(st_kmer, min_size=10, unique_by=lex_rc))
def test_hll_cardinality(kmers):
    """Testing HyperLogLog cardinality estimation.

    len(set) serves as an oracle for HyperLogLog.
    """

    # We use the 'unique_by' argument to guarantee we have only one of
    # {kmer, revcomp(kmer)} in the set: We would overestimate otherwise.

    oracle = set()
    hll = khmer.HLLCounter(0.01, KSIZE)

    for kmer in kmers:
        oracle.update([kmer])
        hll.consume_string(kmer)

    if len(oracle) < 100:
        # This is tricky: if the cardinality is smaller than 100,
        # HyperLogLog might over or underestimate by 1 or 2,
        # which is greater than the default error rate (1%).
        # So, we just check to see if the absolute difference is small.
        assert abs(len(oracle) - len(hll)) <= 2
    else:
        # HyperLogLog error is calculated based on an uniform hash function,
        # but every hash function has some bias, even if small.
        # We set the error rate previously to 1%,
        # but we check for 2% here.
        error = round(abs(len(hll) - len(oracle)) / len(oracle), 2)
        assert error <= 0.02
