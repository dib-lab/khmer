# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2017, The Regents of the University of California.
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
# pylint: disable=missing-docstring,protected-access,no-member,invalid-name

import screed
import khmer
from . import khmer_tst_utils as utils
import pytest


@pytest.mark.parametrize('ksize,memory,epsilon,numbands', [
    (21, 5e6, 1, 2),
    (21, 5e6, 1, 4),
    (21, 5e6, 1, 8),
    (21, 5e6, 1, 16),
])
def test_banding_in_memory(ksize, memory, epsilon, numbands):
    """
    Test accuracy of banding functionally.

    Tests whether k-mer counts loaded into separate counttables in bands gives
    reasonable behavior compared to k-mer counts computed in the normal
    fashion.
    """
    infile = utils.get_test_data('banding-reads.fq')

    ct_normal = khmer.Counttable(ksize, memory / 4, 4)
    ct_normal.consume_seqfile(infile)

    ct_banded = list()
    for band in range(numbands):
        ct = khmer.Counttable(ksize, memory / 4 / numbands, 4)
        ct.consume_seqfile_banding(infile, numbands, band)
        ct_banded.append(ct)

    for n, record in enumerate(screed.open(infile)):
        if not (n > 0 and n % 100 == 0):
            continue
        for kmer in ct_normal.get_kmers(record.sequence):
            abund_normal = ct_normal.get(kmer)
            abunds_banded = [ct.get(kmer) for ct in ct_banded]
            # Ideally, we'd like to enforce that these are equal. Practically,
            # with false positives, we have to allow for a small difference.
            assert abs(sum(abunds_banded) - abund_normal) <= epsilon

            nonzeros = [a for a in abunds_banded if a > 0]
            # False positives shouldn't be appearing in multiple bands
            assert len(nonzeros) <= 2
            # False positives shouldn't have high abundance
            if len(nonzeros) > 1:
                assert min(nonzeros) == 1


@pytest.mark.parametrize('ksize,memory,numbands', [
    (21, 5e6, 3),
    (21, 5e6, 11),
    (21, 5e6, 23),
    (21, 5e6, 29),
])
def test_banding_to_disk(ksize, memory, numbands):
    """
    Test accuracy of banding in terms of the data structure contents.

    Stronger than the functional in-memory test, this function tests whether
    a computing k-mer abundances in banding mode produces the same data
    structure as counting k-mer abundances in the normal fashion.
    """
    infile = utils.get_test_data('banding-reads.fq')
    path1 = utils.get_temp_filename('normal.ct')
    path2 = utils.get_temp_filename('banding.ct')

    ct = khmer.Counttable(ksize, memory / 4, 4)
    ct.consume_seqfile(infile)
    ct.save(path1)
    fpr = khmer.calc_expected_collisions(ct)
    print('FPR', fpr)

    ct = khmer.Counttable(ksize, memory / 4, 4)
    for band in range(numbands):
        ct.consume_seqfile_banding(infile, numbands, band)
    ct.save(path2)
    fpr = khmer.calc_expected_collisions(ct)
    print('FPR', fpr)

    with open(path1, 'rb') as f1, open(path2, 'rb') as f2:
        assert f1.read() == f2.read()


@pytest.mark.parametrize('sketchclass', [
    (khmer.Nodetable),
    (khmer.Counttable),
])
def test_banding_bad_params(sketchclass):
    sketch = sketchclass(31, 1e5, 4)

    # Fails because 13 > 8
    with pytest.raises(ValueError) as ve:
        infile = utils.get_test_data('bogus.fa')
        _ = sketch.consume_seqfile_banding(infile, 8, 13)
    assert "'band' must be in the interval [0, 'num_bands')" in str(ve)

    # Fails because file does not exist
    with pytest.raises(OSError) as ose:
        nreads, kmersconsumed = \
            sketch.consume_seqfile_banding('file-no-exist.fa', 16, 3)
    assert 'does not exist' in str(ose)


@pytest.mark.parametrize('sketchclass,num_batches,batch', [
    (khmer.Nodetable, 8, 3),
    (khmer.Counttable, 8, 3),
])
def test_banding(sketchclass, num_batches, batch):
    sketch = sketchclass(31, 1e5, 4)
    infile = utils.get_test_data('bogus.fa')
    nreads, kmersconsumed = \
        sketch.consume_seqfile_banding(infile, num_batches, batch)
    assert nreads == 1
    assert kmersconsumed == 3
    assert sketch.get('CGGCTATTATCTGAGCTCAAGACTAATACGC') == 1
    assert sketch.get('TATTATCTGAGCTCAAGACTAATACGCGCTG') == 1
    assert sketch.get('TGAGCTCAAGACTAATACGCGCTGGCCACTG') == 1
    assert sketch.get('GTACGGCTATTATCTGAGCTCAAGACTAATA') == 0
    assert sketch.get('TTATCTGAGCTCAAGACTAATACGCGCTGGC') == 0
    assert sketch.get('GCTCAAGACTAATACGCGCTGGCCACTGGTA') == 0
