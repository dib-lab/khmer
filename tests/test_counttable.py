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
#     * Neither the name of the University of California nor the names
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

from __future__ import print_function
from __future__ import absolute_import

import khmer

import pytest
from . import khmer_tst_utils as utils


def test_get_kmer_hashes():
    s = "ATGGATATGGAGGACAAGTATATGGAGGACAAGTATATGGAGGACAAGTAT"
    a = khmer.Counttable(33, 1e6, 3)
    assert a.get_kmer_hashes(s[:33]) == [4743239192574154715]
    assert a.get_kmer_hashes(s[:34]) == [4743239192574154715,
                                         2122462908541313313]

    assert a.get_kmer_hashes(s[0:33]) == [4743239192574154715]
    assert a.get_kmer_hashes(s[1:34]) == [2122462908541313313]


@pytest.mark.parametrize('kmer', [
    ('GATTACA' * 3),
    ('ATG' * 7),
    ('AGGACAAGTATATGGAGGACA'),
])
def test_kmer_revcom_hash(kmer):
    a = khmer.Counttable(21, 1e4, 3)
    assert a.hash(kmer) == a.hash(khmer.reverse_complement(kmer))


@pytest.mark.parametrize('ksize,sketch_allocator', [
    (21, khmer.Nodetable),
    (21, khmer.Counttable),
    (21, khmer.SmallCounttable),
    (49, khmer.Nodetable),
    (49, khmer.Counttable),
    (49, khmer.SmallCounttable),
])
def test_reverse_hash(ksize, sketch_allocator):
    multiplier = int(ksize / len('GATTACA'))
    kmer = 'GATTACA' * multiplier

    sketch = sketch_allocator(ksize, 1e4, 4)
    kmer_hash = sketch.hash(kmer)
    with pytest.raises(ValueError) as ve:
        _ = sketch.reverse_hash(kmer_hash)
    assert 'not implemented' in str(ve)
