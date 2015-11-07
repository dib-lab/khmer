#! /usr/bin/env python

from __future__ import division, unicode_literals

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
