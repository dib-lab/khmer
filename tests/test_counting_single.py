from __future__ import print_function
from __future__ import absolute_import
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#

# pylint: disable=C0111,C0103

import khmer
from . import khmer_tst_utils as utils
from nose.plugins.attrib import attr

MAX_COUNT = 255


def test_no_collision():
    kh = khmer.new_hashtable(4, 4)

    kh.count('AAAA')
    assert kh.get('AAAA') == 1

    kh.count('TTTT')                    # reverse complement
    assert kh.get('TTTT') == 2


@attr('linux')
def test_toobig():
    try:
        ct = khmer.new_hashtable(4, 1000000000000)
        assert 0, "this should fail"
    except MemoryError as err:
        print(str(err))


def test_collision():
    kh = khmer.new_hashtable(4, 4)

    kh.count('AAAA')
    assert kh.get('AAAA') == 1

    kh.count('TTTT')
    assert kh.get('TTTT') == 2


def test_badcount():
    countingtable = khmer.new_hashtable(4, 4)
    try:
        countingtable.count()
        assert 0, "count should require one argument"
    except TypeError as err:
        print(str(err))
    try:
        countingtable.count('ABCDE')
        assert 0, "count should require k-mer size to be equal"
    except ValueError as err:
        print(str(err))


def test_hashtable_n_entries():
    countingtable = khmer.new_hashtable(4, 4)
    try:
        countingtable.n_entries("nope")
        assert 0, "n_entries should accept no arguments"
    except TypeError as err:
        print(str(err))


def test_complete_no_collision():
    kh = khmer.new_hashtable(4, 4 ** 2)

    for i in range(0, kh.n_entries()):
        s = khmer.reverse_hash(i, 4)
        kh.count(s)

    n_palindromes = 0
    n_rc_filled = 0
    n_fwd_filled = 0

    for i in range(0, kh.n_entries()):
        s = khmer.reverse_hash(i, 4)
        if kh.get(s):                   # string hashing is rc aware
            n_rc_filled += 1
        if kh.get(s) == 1:              # palindromes are singular
            n_palindromes += 1
        if kh.get(i):                   # int hashing is not rc aware
            n_fwd_filled += 1

    assert n_rc_filled == kh.n_entries(), n_rc_filled
    assert n_palindromes == 16, n_palindromes  # @CTB check this
    assert n_fwd_filled == kh.n_entries() // 2 + n_palindromes // 2, \
        n_fwd_filled


def test_complete_2_collision():
    kh = khmer.new_hashtable(4, 4)

    for i in range(0, kh.n_entries()):
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
    # @CTB assert n_fwd_filled == 100 # kt.n_entries() / 2, n_fwd_filled


def test_complete_4_collision():
    kh = khmer.new_hashtable(4, 2)

    for i in range(0, kh.n_entries()):
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
    # @CTB assert n_fwd_filled == kt.n_entries() / 2, n_fwd_filled


def test_maxcount():
    # hashtable should saturate at some point so as not to overflow counter
    kh = khmer.new_hashtable(4, 4)

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
    kh = khmer.new_hashtable(4, 4)
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
    kh = khmer.new_hashtable(4, 4)

    s = "TTTT"
    s_rc = "AAAA"

    kh.consume(s)
    n = kh.get(s_rc)
    assert n == 1


def test_maxcount_consume():
    # hashtable should saturate at some point so as not to overflow counter
    kh = khmer.new_hashtable(4, 4)

    s = "A" * 10000
    kh.consume(s)

    c = kh.get('AAAA')
    assert c == MAX_COUNT, c    # this will depend on HashcountType...


def test_maxcount_consume_with_bigcount():
    # use the bigcount hack to avoid saturating the hashtable.
    kh = khmer.new_hashtable(4, 4)
    kh.set_use_bigcount(True)

    s = "A" * 10000
    kh.consume(s)

    c = kh.get('AAAA')
    assert c == 10000 - 3, c


def test_get_mincount():
    kh = khmer.new_hashtable(4, 4)

    s = "AAAAACGT"
    kh.consume(s)

    x = kh.get_min_count(s)
    assert x == 1

    kh.consume(s)
    x = kh.get_min_count(s)
    assert x == 2


def test_get_maxcount():
    kh = khmer.new_hashtable(4, 4)

    s = "AAAAACGT"
    kh.consume(s)

    x = kh.get_max_count(s)
    assert x == 2

    kh.consume(s)
    x = kh.get_max_count(s)
    assert x == 4


def test_get_maxcount_rc():
    kh = khmer.new_hashtable(4, 4)

    s = "AAAAACGT"
    src = "ACGTTTTT"
    kh.consume(s)

    x = kh.get_max_count(s)
    assert x == 2

    kh.consume(src)
    x = kh.get_max_count(s)
    assert x == 4


def test_get_mincount_rc():
    kh = khmer.new_hashtable(4, 4)

    s = "AAAAACGT"
    src = "ACGTTTTT"

    kh.consume(s)
    x = kh.get_min_count(s)
    assert x == 1

    kh.consume(src)
    x = kh.get_min_count(s)
    assert x == 2


def test_badget():
    kh = khmer.new_hashtable(6, 4 ** 10)

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
    kh = khmer.new_hashtable(25, 4)
    fullstr = "GTATGCCAGCTCCAACTGGGCCGGTACGAGCAGGCCATTGCCTCTTGCCGCGATGCGTCGGCG"
    substr = "ATGCCAGCTCCAACTGGGCCGGTACGAGCAGGCCATTGCCTCTTGC"

    kh.consume(fullstr)
    assert 0 < kh.get_min_count(substr), kh.get_min_count(substr)


def test_64bitshift_2():
    kh = khmer.new_hashtable(25, 4)
    fullstr = "GTATGCCAGCTCCAACTGGGCCGGTACGAGCAGGCCATTGCCTCTTGCCGCGATGCGTCGGCG"

    kh.consume(fullstr)
    for i in range(len(fullstr) - 25 + 1):
        substr = fullstr[i:i + 25]
        assert kh.get(substr) > 0


def test_very_short_read():
    short_filename = utils.get_test_data('test-short.fa')
    kh = khmer.new_hashtable(9, 4)
    n_reads, n_kmers = kh.consume_fasta(short_filename)
    assert n_reads == 1, n_reads
    assert n_kmers == 0, n_kmers

    kh = khmer.new_hashtable(8, 4)
    n_reads, n_kmers = kh.consume_fasta(short_filename)
    assert n_reads == 1, n_reads
    assert n_kmers == 1, n_kmers


class Test_ConsumeString(object):

    def setup(self):
        self.kh = khmer.new_hashtable(4, 4 ** 4)

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

    def test_abundance_by_pos(self):
        kh = self.kh

        for _ in range(0, 300):
            kh.count('ATCG')

        for _ in range(0, 10):
            kh.count('ATGG')

        short_filename = utils.get_test_data('test-short.fa')
        dist = kh.fasta_count_kmers_by_position(short_filename, 6, 10)
        assert dist[4] == 1
        assert sum(dist) == 1

        dist = kh.fasta_count_kmers_by_position(short_filename, 6, MAX_COUNT)
        assert dist[0] == 1, dist[0]
        assert dist[2] == 1
        assert sum(dist) == 2

    def test_abundance_by_pos_bigcount(self):
        kh = self.kh
        kh.set_use_bigcount(True)       # count past MAX_COUNT

        for _ in range(0, 300):
            kh.count('ATCG')

        for _ in range(0, 10):
            kh.count('ATGG')

        short_filename = utils.get_test_data('test-short.fa')
        dist = kh.fasta_count_kmers_by_position(short_filename, 6, 10)
        assert dist[4] == 1
        assert sum(dist) == 1

        dist = kh.fasta_count_kmers_by_position(short_filename, 6, 300)
        assert dist[0] == 1, dist[0]
        assert dist[2] == 1
        assert sum(dist) == 2

    def test_n_occupied_args(self):
        assert self.kh.n_occupied() == 0
        self.kh.consume('AAAA')
        assert self.kh.n_occupied(0, 1) == 1
        assert self.kh.n_occupied(1, 4 ** 4) == 0

        hashvalue = khmer.forward_hash('AACT', 4)
        self.kh.consume('AACT')
        assert self.kh.n_occupied(0, hashvalue + 1) == 2
        assert self.kh.n_occupied(hashvalue + 1, 4 ** 4) == 0
        assert self.kh.n_occupied(hashvalue, hashvalue + 1) == 1

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
        self.kh = khmer.new_hashtable(4, 4)
        A_filename = utils.get_test_data('all-A.fa')
        self.kh.consume_fasta(A_filename)

    def test_count_A(self):
        A_filename = utils.get_test_data('all-A.fa')

        tracking = khmer.new_hashbits(4, 4, 1)
        dist = self.kh.abundance_distribution(A_filename, tracking)

        assert sum(dist) == 1
        assert dist[10] == 1
