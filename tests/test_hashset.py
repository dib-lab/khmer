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
# pylint: disable=C0111,C0103
"""
Test code for HashSet objects.
"""

import khmer
from . import khmer_tst_utils as utils


def test_bad_construct():
    try:
        hs = khmer.HashSet()
        assert 0, "HashSet constructor should fail w/o argument"
    except TypeError:
        pass

    try:
        hs = khmer.HashSet(5, [{}])
        assert 0, "HashSet constructor should fail w/o list of k-mers"
    except TypeError:
        pass


def test_iter_single():
    hs = khmer.HashSet(5, [6])
    for k in hs:
        assert k == 6
        print(k)


def test_iter_double():
    x = [6, 9, 20]
    hs = khmer.HashSet(5, x)
    for i, k in enumerate(hs):
        assert k == x[i], (k, x[i])


def test_iter_single():
    hs = khmer.HashSet(5, [6])
    k = iter(hs)
    k2 = iter(k)
    assert k == k2


def test_add():
    hs = khmer.HashSet(5)
    hs.add(7)
    hs.add(4)

    assert list(sorted(hs)) == [4, 7]


def test_update():
    hs = khmer.HashSet(5)
    x = [5, 10, 15, 2**35]
    hs.update(x)

    assert list(sorted(hs)) == [5, 10, 15, 2**35]


def test_update_bad():
    hs = khmer.HashSet(5)
    x = [5, 10, 15, 2**35, {}]
    try:
        hs.update(x)
        assert 0, "cannot add dict to a HashSet"
    except TypeError:
        pass


def test_remove():
    hs = khmer.HashSet(5, [8, 10])
    assert len(hs) == 2
    hs.remove(8)
    assert len(hs) == 1
    assert list(hs) == [10]


def test_remove_2():
    hs = khmer.HashSet(5, [8, 10])
    assert len(hs) == 2
    try:
        hs.remove(15)
        assert 0, "hs.remove should raise an Exception"
    except ValueError:
        pass
    assert len(hs) == 2
    assert list(sorted(hs)) == [8, 10]


def test_contains_1():
    hs = khmer.HashSet(5, [8, 10])
    assert 8 in hs
    assert 10 in hs
    assert 2**35 not in hs


def test_contains_2():
    hs = khmer.HashSet(5, [8, 10])
    assert khmer.reverse_hash(8, 5) in hs
    assert khmer.reverse_hash(10, 5) in hs
    assert khmer.reverse_hash(2**35, 5) not in hs


def test_concat_1():
    hs = khmer.HashSet(5, [10, 12])
    hs2 = khmer.HashSet(5, [10, 13])

    hs3 = hs + hs2
    assert list(sorted(hs3)) == [10, 12, 13]


def test_concat_2():
    hs = khmer.HashSet(5, [10, 12])
    hs2 = khmer.HashSet(5, [10, 13])

    hs += hs2
    assert list(sorted(hs)) == [10, 12, 13]


def test_concat_1_fail():
    hs = khmer.HashSet(5, [10, 12])
    hs2 = khmer.HashSet(4, [10, 13])

    try:
        hs3 = hs + hs2
        assert 0, "concat should fail - different ksize"
    except ValueError:
        pass


def test_concat_2_fail():
    hs = khmer.HashSet(5, [10, 12])
    hs2 = khmer.HashSet(4, [10, 13])

    try:
        hs += hs2
        assert 0, "inplace concat should fail - different ksize"
    except ValueError:
        pass
