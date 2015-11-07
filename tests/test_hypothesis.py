#! /usr/bin/env python

from __future__ import division, unicode_literals

from collections import Counter

from hypothesis import given, assume, strategies as st
from nose.plugins.attrib import attr

import khmer
from khmer import reverse_hash, forward_hash, forward_hash_no_rc


TRANSLATE = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}

KSIZE = 13
st_kmer = st.text("ACGT", min_size=KSIZE, max_size=KSIZE)


def lex_rc(kmer):
    return reverse_hash(forward_hash(kmer, len(kmer)), len(kmer))


@attr('hypothesis')
@given(st_kmer)
def test_forward_hash(kmer):
    ksize = len(kmer)

    rh = reverse_hash(forward_hash(kmer, ksize), ksize)
    assert rh == kmer or rh == "".join(TRANSLATE[c] for c in kmer[::-1])


@attr('hypothesis')
@given(st_kmer)
def test_forward_hash_no_rc(kmer):
    ksize = len(kmer)

    rh = reverse_hash(forward_hash_no_rc(kmer, ksize), ksize)
    assert rh == kmer


@attr('hypothesis')
@given(st.lists(st_kmer, min_size=1))
def test_countgraph_undercounting(kmers):
    oracle = Counter()
    countgraph = khmer.Countgraph(KSIZE, 4 ** 4, 4)

    for kmer in kmers:
        oracle.update([kmer])
        countgraph.count(kmer)

    for kmer in oracle:
        assert countgraph.get(kmer) >= min(oracle[kmer], 255)


@attr('hypothesis')
@given(st.sets(st_kmer, min_size=1))
def test_nodegraph_presence(kmers):
    oracle = set()
    nodegraph = khmer.Nodegraph(KSIZE, 4 ** 4, 4)

    for kmer in kmers:
        oracle.update([kmer])
        nodegraph.count(kmer)

    for kmer in oracle:
        assert nodegraph.get(kmer) == 1
