#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import sys
import khmer
import string
from array import array

__complementTranslation = string.maketrans('ACTG', 'TGAC')


def complement(s):
    """
    Return complement of 's'.
    """
    c = string.translate(s, __complementTranslation)
    return c

#
# reverse
#


def reverse(s):
    """
    Return reverse of 's'.
    """
    r = array('c', s)
    r.reverse()
    r = string.join(r, '')

    return r


def rc(s):
    return reverse(complement(s))


def test_rc():
    assert rc("ACTG") == "CAGT"

# build a k-table of length L, and test it.

L = 4


class Test_KTable:

    def setup(self):
        # make a new ktable.
        self.kt = khmer.new_ktable(L)

    def test_basic(self):
        kt = self.kt
        # check to make sure sizes are what we expect.
        assert kt.ksize() == L
        assert kt.max_hash() == 4 ** L - 1
        assert kt.n_entries() == 4 ** L
        assert len(kt) == 4 ** L

    def test_hash(self):
        kt = self.kt

        # make sure forward/reverse hash work minimally.
        s = 'ATCG'
        assert kt.reverse_hash(kt.forward_hash('ATCG')) == s

    def test_populate(self):
        return                          # @CTB
        kt = self.kt

        # make sure hashes populate the table completely, too!

        alphabet = ('A', 'C', 'G', 'T')

        def rN(L, d={}, *args):
            if L == 0:
                d["".join(args)] = 1
                return

            for letter in alphabet:
                rN(L - 1, d, letter, *args)

            return d.keys()

        # generate all L-mers & make sure they all map differently:
        all_lmers = rN(L)

        occupy_l = []
        for i in all_lmers:
            occupy_l.append(kt.forward_hash(i))

        occupy_l.sort()
        assert occupy_l == range(0, kt.n_entries())

        # check to make sure that fwd --> rev --> fwd works.

        for i in range(0, kt.n_entries()):
            assert kt.forward_hash(kt.reverse_hash(i)) == i

    def test_consume(self):
        return                          # @CTB
        kt = self.kt

        # consume a test string, and verify that consume works.
        s = "ATGAGAGACACAGGGAGAGACCCAATTAGAGAATTGGACC"
        kt.consume(s)

        kt2 = khmer.new_ktable(L)

        for start in range(0, len(s) - L + 1):
            word = s[start:start + L]

            kt2.count(word)

        for i in range(0, kt.n_entries()):
            n = kt.get(i)                       # test 'consume_str' numbers
            n3 = kt2.get(i)                     # and 'count' count.
            assert n == n3

        for i in range(0, kt.n_entries()):
            kt.set(i, 1)

        for i in range(0, kt.n_entries()):
            assert(kt.get(i) == 1)

    def test_operator_in(self):
        kt = self.kt

        s = "ATGAGAGACACAGGGAGAGACCCAATTAGAGAATTGGACC"
        kt.consume(s)

        assert "CCCAA" in kt
        assert "GGGGG" not in kt

    def test_intersection(self):
        kt = self.kt

        # intersection
        for i in range(0, 4 ** L / 4):
            kt.set(i * 4, 1)

        kt2 = khmer.new_ktable(L)
        for i in range(0, 4 ** L / 5):
            kt2.set(i * 5, 1)

        kt3 = kt.intersect(kt2)

        assert kt3.get(20) == 2
        for i in range(0, 4 ** L):
            if kt3.get(i):
                assert i % 4 == 0
                assert i % 5 == 0

    def test_update(self):
        kt = self.kt

        # intersection
        for i in range(0, 4 ** L / 4):
            kt.set(i * 4, 1)

        kt2 = khmer.new_ktable(L)
        for i in range(0, 4 ** L / 5):
            kt2.set(i * 5, 1)

        kt.update(kt2)
        for i in range(0, 4 ** L):
            if kt.get(i):
                assert i % 4 == 0 or i % 5 == 0

    def test_clear(self):
        kt = self.kt

        # test clear()
        for i in range(0, 4 ** L / 4):
            kt.set(i * 4, 1)

        kt.clear()
        for i in range(0, 4 ** L):
            assert kt.get(i) == 0
            assert kt[i] == 0
