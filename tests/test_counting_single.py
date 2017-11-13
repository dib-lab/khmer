# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2010-2015, Michigan State University.
# Copyright (C) 2015-2016, The Regents of the University of California.
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
# pylint: disable=C0111,C0103,missing-docstring,no-member,protected-access


import khmer

import pytest
from . import khmer_tst_utils as utils


MAX_COUNT = 255


@pytest.mark.huge
def test_toobig():
    try:
        khmer.Countgraph(4, 1000000000000, 1)
        assert 0, "this should fail"
    except MemoryError as err:
        print(str(err))


def test_collision():
    kh = khmer.Countgraph(4, 100, 1)

    kh.count('AAAA')
    assert kh.get('AAAA') == 1

    kh.count('TTTT')
    assert kh.get('TTTT') == 2


def test_badcount():
    countgraph = khmer.Countgraph(4, 1, 1)
    try:
        countgraph.count()
        assert 0, "count should require one argument"
    except TypeError as err:
        print(str(err))
    try:
        countgraph.count('ABCDE')
        assert 0, "count should require k-mer size to be equal"
    except ValueError as err:
        print(str(err))


def test_complete_no_collision():
    kh = khmer.Countgraph(4, 1, 1, primes=[4 ** 4])

    n_entries = kh.hashsizes()[0]

    for i in range(0, n_entries):
        s = khmer.reverse_hash(i, 4)
        kh.count(s)

    n_palindromes = 0
    n_rc_filled = 0
    n_fwd_filled = 0

    for i in range(0, n_entries):
        s = khmer.reverse_hash(i, 4)
        if kh.get(s):                   # string hashing is rc aware
            n_rc_filled += 1
        if kh.get(s) == 1:              # palindromes are singular
            n_palindromes += 1
        if kh.get(i):                   # int hashing is not rc aware
            n_fwd_filled += 1

    assert n_rc_filled == n_entries, n_rc_filled
    assert n_palindromes == 16, n_palindromes
    assert n_fwd_filled == n_entries // 2 + n_palindromes // 2, \
        (n_fwd_filled, n_entries // 2 + n_palindromes // 2)


def test_complete_2_collision():
    kh = khmer.Countgraph(4, 7, 1)

    n_entries = kh.hashsizes()[0]
    for i in range(0, n_entries):
        s = khmer.reverse_hash(i, 4)
        kh.count(s)

    n_rc_filled = 0
    #  n_fwd_filled = 0

    for i in range(0, 128):
        s = khmer.reverse_hash(i, 4)
        if kh.get(s):                   # string hashing is rc aware
            n_rc_filled += 1
    # if kh.get(i):                   # int hashing is not rc aware
    #        n_fwd_filled += 1

    assert n_rc_filled == 128, n_rc_filled


def test_complete_4_collision():
    kh = khmer.Countgraph(4, 5, 1)

    n_entries = kh.hashsizes()[0]

    for i in range(0, n_entries):
        s = khmer.reverse_hash(i, 4)
        kh.count(s)

    n_rc_filled = 0
    #  n_fwd_filled = 0

    for i in range(0, 64):
        s = khmer.reverse_hash(i, 4)
        if kh.get(s):                   # string hashing is rc aware
            n_rc_filled += 1
    # if kh.get(i):                   # int hashing is not rc aware
    #       n_fwd_filled += 1

    assert n_rc_filled == 64, n_rc_filled


def test_maxcount():
    # hashtable should saturate at some point so as not to overflow counter
    kh = khmer.Countgraph(4, 100, 1)

    last_count = None
    for _ in range(0, 10000):
        kh.count('AAAA')
        c = kh.get('AAAA')

        print(last_count, c)
        if c == last_count:
            break
        last_count = c

    assert c != 10000, "should not be able to count to 10000"
    assert c == MAX_COUNT       # this will depend on HashcountType...


def test_maxcount_with_bigcount():
    # hashtable should not saturate, if use_bigcount is set.
    kh = khmer.Countgraph(4, 100, 1)
    kh.set_use_bigcount(True)

    last_count = None
    for _ in range(0, 10000):
        kh.count('AAAA')
        c = kh.get('AAAA')

        print(last_count, c)
        if c == last_count:
            break
        last_count = c

    assert c == 10000, "should be able to count to 10000"
    assert c != MAX_COUNT


def test_consume_uniqify_first():
    kh = khmer.Countgraph(4, 100, 1)

    s = "TTTT"
    s_rc = "AAAA"

    kh.consume(s)
    n = kh.get(s_rc)
    assert n == 1


def test_maxcount_consume():
    # hashtable should saturate at some point so as not to overflow counter
    kh = khmer.Countgraph(4, 100, 1)

    s = "A" * 10000
    kh.consume(s)

    c = kh.get('AAAA')
    assert c == MAX_COUNT, c    # this will depend on HashcountType...


def test_maxcount_consume_with_bigcount():
    # use the bigcount hack to avoid saturating the hashtable.
    kh = khmer.Countgraph(4, 100, 1)
    kh.set_use_bigcount(True)

    s = "A" * 10000
    kh.consume(s)

    c = kh.get('AAAA')
    assert c == 10000 - 3, c


def test_get_mincount():
    kh = khmer.Countgraph(4, 100, 1)

    s = "AAAAACGT"
    kh.consume(s)

    x = kh.get_min_count(s)
    assert x == 1, x

    kh.consume(s)
    x = kh.get_min_count(s)
    assert x == 2, x


def test_get_maxcount():
    kh = khmer.Countgraph(4, 9, 1)

    s = "AAAAACGT"
    kh.consume(s)

    x = kh.get_max_count(s)
    assert x == 2

    kh.consume(s)
    x = kh.get_max_count(s)
    assert x == 4


def test_get_maxcount_rc():
    kh = khmer.Countgraph(4, 9, 1)

    s = "AAAAACGT"
    src = "ACGTTTTT"
    kh.consume(s)

    x = kh.get_max_count(s)
    assert x == 2, x

    kh.consume(src)
    x = kh.get_max_count(s)
    assert x == 4, x


def test_get_mincount_rc():
    kh = khmer.Countgraph(4, 100, 1)

    s = "AAAAACGT"
    src = "ACGTTTTT"

    kh.consume(s)
    x = kh.get_min_count(s)
    assert x == 1, x

    kh.consume(src)
    x = kh.get_min_count(s)
    assert x == 2


def test_badget():
    kh = khmer.Countgraph(6, 4 ** 10, 1)

    DNA = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAG"

    kh.consume(DNA)

    assert kh.get("AGCTTT") == 1

    assert kh.get("GATGAG") == 0

    try:
        kh.get("AGCTT")
        assert 0, "this should fail"
    except ValueError as err:
        print(str(err))


def test_64bitshift():
    kh = khmer.Countgraph(25, 4, 1)
    fullstr = "GTATGCCAGCTCCAACTGGGCCGGTACGAGCAGGCCATTGCCTCTTGCCGCGATGCGTCGGCG"
    substr = "ATGCCAGCTCCAACTGGGCCGGTACGAGCAGGCCATTGCCTCTTGC"

    kh.consume(fullstr)
    assert 0 < kh.get_min_count(substr), kh.get_min_count(substr)


def test_64bitshift_2():
    kh = khmer.Countgraph(25, 4, 1)
    fullstr = "GTATGCCAGCTCCAACTGGGCCGGTACGAGCAGGCCATTGCCTCTTGCCGCGATGCGTCGGCG"

    kh.consume(fullstr)
    for i in range(len(fullstr) - 25 + 1):
        substr = fullstr[i:i + 25]
        assert kh.get(substr) > 0


def test_very_short_read():
    short_filename = utils.get_test_data('test-short.fa')
    kh = khmer.Countgraph(9, 4, 1)
    n_reads, n_kmers = kh.consume_seqfile(short_filename)
    assert n_reads == 1, n_reads
    assert n_kmers == 0, n_kmers

    kh = khmer.Countgraph(8, 4, 1)
    n_reads, n_kmers = kh.consume_seqfile(short_filename)
    assert n_reads == 1, n_reads
    assert n_kmers == 1, n_kmers


class Test_ConsumeString(object):

    def setup(self):
        self.kh = khmer.Countgraph(4, 1, 1, primes=[4 ** 4])

    def test_n_occupied(self):
        assert self.kh.n_occupied() == 0
        self.kh.consume('AAAA')
        assert self.kh.n_occupied() == 1
        self.kh.consume('AACT')
        assert self.kh.n_occupied() == 2
        try:
            self.kh.n_occupied("MU", 1, 3)
            assert 0, "n_occupied shouldn't accept three arguments"
        except TypeError as err:
            print(str(err))

    def test_simple(self):
        n = self.kh.consume('AAAA')
        assert n == 1
        assert self.kh.get(0) == 1

    def test_simple_2(self):
        n = self.kh.consume('AAAAA')
        assert n == 2
        assert self.kh.get(0) == 2

    def test_simple_rc(self):
        n = self.kh.consume('TTTTT')
        assert n == 2
        assert self.kh.get(0) == 2

    def test_min_count(self):
        self.kh.consume('AAAA')

        count = self.kh.get_min_count('AAAA')
        assert count == 1

    def test_max_count(self):
        self.kh.consume('AAAA')

        count = self.kh.get_max_count('AAAA')
        assert count == 1


class Test_AbundanceDistribution(object):

    def setup(self):
        self.kh = khmer.Countgraph(4, 100, 1)
        A_filename = utils.get_test_data('all-A.fa')
        self.kh.consume_seqfile(A_filename)

    def test_count_A(self):
        A_filename = utils.get_test_data('all-A.fa')

        tracking = khmer.Nodegraph(4, 7, 1)
        dist = self.kh.abundance_distribution(A_filename, tracking)

        assert sum(dist) == 1
        assert dist[10] == 1
