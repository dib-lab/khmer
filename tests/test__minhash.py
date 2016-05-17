# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2016, The Regents of the University of California.
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
# Contact: titus@idyll.org
# pylint: disable=missing-docstring,protected-access

from __future__ import print_function
from __future__ import absolute_import, unicode_literals

from khmer import MinHash
import screed

# add:
# * get default params from Python
# * keyword args for minhash constructor
# * trap error from handing protein/non-DNA to a DNA MH
# * fail on untagged/unloaded countgraph
# * nan on empty minhash
# * define equals
# * check prime, ksize, etc before comparing and merging

def test_default_params():
    # verify that MHs have these default parameters.
    mh = MinHash(100, 32)
    mh2 = MinHash(100, 32, 9999999967, False)  # should not be changed!

    mh.add_hash(9999999968)
    mh2.add_hash(9999999968)
    assert mh.get_mins() == mh2.get_mins()
    assert mh.get_mins()[0] == 1

def test_basic_dna():
    # verify that MHs of size 1 stay size 1, & act properly as bottom sketches.
    mh = MinHash(1, 4)
    mh.add_sequence('ATGC')
    a = mh.get_mins()

    mh.add_sequence('GCAT')             # this will not get added; hash > ATGC
    b = mh.get_mins()

    print(a, b)
    assert a == b
    assert len(b) == 1


def test_protein():
    # verify that we can hash protein/aa sequences
    mh = MinHash(1, 4, 9999999967, True)
    mh.add_protein('AGYYG')

    assert len(mh.get_mins()) == 1


def test_compare_1():
    a = MinHash(20, 10)
    b = MinHash(20, 10)

    a.add_sequence('TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA')
    b.add_sequence('TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA')

    assert a.compare(b) == 1.0
    assert b.compare(b) == 1.0
    assert b.compare(a) == 1.0
    assert a.compare(a) == 1.0

    # add same sequence again
    b.add_sequence('TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA')
    assert a.compare(b) == 1.0
    assert b.compare(b) == 1.0
    assert b.compare(a) == 1.0
    assert a.compare(a) == 1.0

    
    b.add_sequence('GATTGGTGCACACTTAACTGGGTGCCGCGCTGGTGCTGATCCATGAAGTT')
    x = a.compare(b)
    assert x >= 0.3, x
    
    x = b.compare(a)
    assert x >= 0.3, x
    assert a.compare(a) == 1.0
    assert b.compare(b) == 1.0


def test_mh_copy():
    a = MinHash(20, 10)

    a.add_sequence('TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA')
    b = a.__copy__()
    assert b.compare(a) == 1.0


def test_mh_len():
    a = MinHash(20, 10)

    assert len(a) == 20
    a.add_sequence('TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA')
    assert len(a) == 20


def test_mh_len():
    a = MinHash(20, 10)
    for i in range(0, 40, 2):
        a.add_hash(i)

    assert a.get_mins() == list(range(0, 40, 2))


def test_mh_count_common():
    a = MinHash(20, 10)
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = MinHash(20, 10)
    for i in range(0, 80, 4):
        b.add_hash(i)

    assert a.count_common(b) == 10
    assert b.count_common(a) == 10


def test_mh_count_common_diff_prime():
    a = MinHash(20, 5, 9999999967)
    b = MinHash(20, 5, 9999999968)

    try:
        a.count_common(b)
        assert 0, "count_common should fail with different prime"
    except ValueError:
        pass


def test_mh_count_common_diff_protein():
    a = MinHash(20, 5, 9999999967, False)
    b = MinHash(20, 5, 9999999967, True)

    try:
        a.count_common(b)
        assert 0, "count_common should fail with DNA vs protein"
    except ValueError:
        pass


def test_mh_count_common_diff_ksize():
    a = MinHash(20, 5)
    b = MinHash(20, 6)

    try:
        a.count_common(b)
        assert 0, "count_common should fail with different ksize"
    except ValueError:
        pass


def test_mh_asymmetric():
    a = MinHash(20, 10)
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = MinHash(10, 10)                   # different size: 10
    for i in range(0, 80, 4):
        b.add_hash(i)

    assert a.count_common(b) == 10
    assert b.count_common(a) == 10

    assert a.compare(b) == 0.5
    assert b.compare(a) == 1.0


def test_mh_merge():
    # test merging two identically configured minhashes
    a = MinHash(20, 10)
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = MinHash(20, 10)
    for i in range(0, 80, 4):
        b.add_hash(i)

    c = a.merge(b)
    d = b.merge(a)

    assert len(c) == len(d)
    assert c.get_mins() == d.get_mins()
    assert c.compare(d) == 1.0
    assert d.compare(c) == 1.0




def test_mh_asymmetric_merge():
    # test merging two asymmetric (different size) MHs
    a = MinHash(20, 10)
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = MinHash(10, 10)                   # different size: 10
    for i in range(0, 80, 4):
        b.add_hash(i)

    c = a.merge(b)
    d = b.merge(a)

    assert len(a) == 20
    assert len(b) == 10
    assert len(c) == len(a)
    assert len(d) == len(b)

    assert d.compare(a) == 1.0
    assert c.compare(b) == 0.5


def test_mh_inplace_concat_asymmetric():
    # test merging two asymmetric (different size) MHs
    a = MinHash(20, 10)
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = MinHash(10, 10)                   # different size: 10
    for i in range(0, 80, 4):
        b.add_hash(i)

    c = a.__copy__()
    c += b

    d = b.__copy__()
    d += a

    assert len(a) == 20
    assert len(b) == 10
    assert len(c) == len(a)
    assert len(d) == len(b)

    assert d.compare(a) == 1.0
    assert c.compare(b) == 0.5


def test_mh_inplace_concat():
    # test merging two identically configured minhashes
    a = MinHash(20, 10)
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = MinHash(20, 10)
    for i in range(0, 80, 4):
        b.add_hash(i)

    c = a.__copy__()
    c += b
    d = b.__copy__()
    d += a

    assert len(c) == len(d)
    assert c.get_mins() == d.get_mins()
    assert c.compare(d) == 1.0
    assert d.compare(c) == 1.0

def test_mh_merge_diff_prime():
    a = MinHash(20, 5, 9999999967)
    b = MinHash(20, 5, 9999999968)

    try:
        a.merge(b)
        assert 0, "merge should fail with different prime"
    except ValueError:
        pass


def test_mh_merge_diff_protein():
    a = MinHash(20, 5, 9999999967, False)
    b = MinHash(20, 5, 9999999967, True)

    try:
        a.merge(b)
        assert 0, "merge should fail with DNA vs protein"
    except ValueError:
        pass


def test_mh_merge_diff_ksize():
    a = MinHash(20, 5)
    b = MinHash(20, 6)

    try:
        a.merge(b)
        assert 0, "merge should fail with different ksize"
    except ValueError:
        pass


def test_mh_compare_diff_prime():
    a = MinHash(20, 5, 9999999967)
    b = MinHash(20, 5, 9999999968)

    try:
        a.compare(b)
        assert 0, "compare should fail with different prime"
    except ValueError:
        pass


def test_mh_compare_diff_protein():
    a = MinHash(20, 5, 9999999967, False)
    b = MinHash(20, 5, 9999999967, True)

    try:
        a.compare(b)
        assert 0, "compare should fail with DNA vs protein"
    except ValueError:
        pass


def test_mh_compare_diff_ksize():
    a = MinHash(20, 5)
    b = MinHash(20, 6)

    try:
        a.compare(b)
        assert 0, "compare should fail with different ksize"
    except ValueError:
        pass


def test_mh_concat_diff_prime():
    a = MinHash(20, 5, 9999999967)
    b = MinHash(20, 5, 9999999968)

    try:
        a += b
        assert 0, "concat should fail with different prime"
    except ValueError:
        pass


def test_mh_concat_diff_protein():
    a = MinHash(20, 5, 9999999967, False)
    b = MinHash(20, 5, 9999999967, True)

    try:
        a += b
        assert 0, "concat should fail with DNA vs protein"
    except ValueError:
        pass


def test_mh_concat_diff_ksize():
    a = MinHash(20, 5)
    b = MinHash(20, 6)

    try:
        a += b
        assert 0, "concat should fail with different ksize"
    except ValueError:
        pass
