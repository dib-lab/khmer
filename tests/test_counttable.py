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


import khmer
import pytest
import random
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


def test_consume_with_mask():
    """
    Test bulk loading with a mask

    The top sequence is the mask, the bottom sequence is to be loaded. The
    bottom 3 k-mers are not present in the mask and therefore should be the
    only ones loaded into the counttable.

    TAGATCTGCTTGAAACAAGTGGATTTGAGAAAAA        <--- mask
       ATCTGCTTGAAACAAGTGGATTTGAGAAAAAAGT     <--- sequence
                          |-----------|
                           |-----------|
                            |-----------|
    """
    maskfile = utils.get_test_data('seq-a.fa')
    mask = khmer.Counttable(13, 1e3, 4)
    mask.consume_seqfile(maskfile)

    infile = utils.get_test_data('seq-b.fa')
    ct = khmer.Counttable(13, 1e3, 4)
    nr, nk = ct.consume_seqfile_with_mask(infile, mask)

    assert nr == 1
    assert nk == 3
    assert ct.get('GATTTGAGAAAAA') == 0  # in the mask
    assert ct.get('ATTTGAGAAAAAA') == 1
    assert ct.get('TTTGAGAAAAAAG') == 1
    assert ct.get('TTGAGAAAAAAGT') == 1


def test_consume_banding_with_mask():
    """
    Test bulk loading with a mask *in k-mer banding mode*

    TAGATCTGCTTGAAACAAGTGGATTTGAGAAAAA        <--- mask
       ATCTGCTTGAAACAAGTGGATTTGAGAAAAAAGT     <--- sequence
                          |-----------|
                           |-----------|
                            |-----------|     <--- only k-mer in band 1/4
    """
    maskfile = utils.get_test_data('seq-a.fa')
    mask = khmer.Counttable(13, 1e3, 4)
    mask.consume_seqfile(maskfile)

    infile = utils.get_test_data('seq-b.fa')
    ct = khmer.Counttable(13, 1e3, 4)
    nr, nk = ct.consume_seqfile_banding_with_mask(infile, 4, 1, mask)

    assert nr == 1
    assert nk == 1
    assert ct.get('GATTTGAGAAAAA') == 0  # in the mask
    assert ct.get('ATTTGAGAAAAAA') == 0  # out of band
    assert ct.get('TTTGAGAAAAAAG') == 0  # out of band
    assert ct.get('TTGAGAAAAAAGT') == 1


def test_consume_with_mask_threshold():
    """
    Test bulk loading with a mask and an abundance threshold

    The top sequence is the mask, the bottom sequence is to be loaded. The
    bottom 3 k-mers are not present in the mask and therefore should be the
    only ones loaded into the counttable.

    TAGATCTGCTTGAAACAAGTGGATTTGAGAAAAA        |
    TAGATCTGCTTGAAACAAGTGGATTTGAGAAAAA        |
    TAGATCTGCTTGAAACAAGTGGATTTGAGAAAAA        | <--- mask input
    TAGATCTGCTTGAAACAAGTGGATTTGAGAAAAAAGT     |
    TAGATCTGCTTGAAACAAGTGGATTTGAGAAAAAAGT     |

       ATCTGCTTGAAACAAGTGGATTTGAGAAAAAAGT     <--- sequence
                          |-----------|
                           |-----------|      <--- only these k-mers are
                            |-----------|          abundance <= 3 in the mask
    """
    mask = khmer.Counttable(13, 1e3, 4)
    for _ in range(3):
        mask.consume('TAGATCTGCTTGAAACAAGTGGATTTGAGAAAAA')
    for _ in range(2):
        mask.consume('TAGATCTGCTTGAAACAAGTGGATTTGAGAAAAAAGT')

    infile = utils.get_test_data('seq-b.fa')
    ct = khmer.Counttable(13, 1e3, 4)
    nr, nk = ct.consume_seqfile_with_mask(infile, mask, 3)

    assert nr == 1
    assert nk == 3
    assert ct.get('GATTTGAGAAAAA') == 0  # in the mask
    assert ct.get('ATTTGAGAAAAAA') == 1
    assert ct.get('TTTGAGAAAAAAG') == 1
    assert ct.get('TTGAGAAAAAAGT') == 1


def test_consume_with_mask_complement():
    mask = khmer.Nodetable(13, 1e3, 4)
    mask.consume('TGCTTGAAACAAGTG')

    infile = utils.get_test_data('seq-b.fa')
    ct = khmer.Counttable(13, 1e3, 4)
    nr, nk = ct.consume_seqfile_with_mask(infile, mask, threshold=1,
                                          consume_masked=True)

    assert ct.get_kmer_counts('TGCTTGAAACAAGTG') == [1, 1, 1]
    assert ct.get_kmer_counts('GAAACAAGTGGATTT') == [0, 0, 0]


@pytest.mark.parametrize('sketchtype', [
    (khmer.Nodegraph),
    (khmer.Countgraph),
    (khmer.SmallCountgraph),
    (khmer.Nodetable),
    (khmer.Counttable),
    (khmer.SmallCounttable),
    (khmer.CyclicCounttable),
])
def test_init_with_primes(sketchtype):
    primes = khmer.get_n_primes_near_x(4, random.randint(1000, 2000))
    sketch = sketchtype(31, 1, 1, primes=primes)
    assert sketch.hashsizes() == primes
