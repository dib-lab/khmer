# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2010-2015, Michigan State University.
# Copyright (C) 2015-2016, The Regents of the University of California.
# Copyright (C) 2016, Google, Inc
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

import gzip

import os

import khmer
from khmer import QFCounttable
from . import khmer_tst_utils as utils
from khmer import ReadParser
import screed
import random
import pytest

MAX_COUNT = 255
MAX_BIGCOUNT = 65535

sketchSize=1048576
 

DNA = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"




def teardown():
    utils.cleanup()





@pytest.fixture(params=[khmer.QFCounttable])
def getSketch(request):
    return request.param

    

def test_count_1(getSketch): 
    hi = getSketch(12, sketchSize)
        
    kmer = 'G' * 12
    hashval = hi.hash('G' * 12)
        
    assert hi.get(kmer) == 0
    assert hi.get(hashval) == 0
        
    hi.count(kmer)
    assert hi.get(kmer) == 1
    assert hi.get(hashval) == 1

    hi.count(kmer)
    assert hi.get(kmer) == 2
    assert hi.get(hashval) == 2
        
    kmer = 'G' * 11
        
    with pytest.raises(ValueError):
        hi.hash(kmer)

            
def test_count_2(getSketch):
    hi = getSketch(12, sketchSize)
    print("done")
    kmer = 'G' * 12
    hashval = hi.hash('G' * 12)

    assert hi.get(kmer) == 0
    assert hi.get(hashval) == 0

    hi.count(kmer)
    assert hi.get(kmer) == 1
    assert hi.get(hashval) == 1

    hi.count(hashval)                     # count hashes same as strings
    assert hi.get(kmer) == 2
    assert hi.get(hashval) == 2



def test_consume_with_mask(getSketch):
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
    mask = getSketch(13, sketchSize)
    mask.consume_seqfile(maskfile)

    infile = utils.get_test_data('seq-b.fa')
    ct =getSketch(12, sketchSize)
    nr, nk = ct.consume_seqfile_with_mask(infile, mask)

    assert nr == 1
    assert nk == 3
    assert ct.get('GATTTGAGAAAAA') == 0  # in the mask
    assert ct.get('ATTTGAGAAAAAA') == 1
    assert ct.get('TTTGAGAAAAAAG') == 1
    assert ct.get('TTGAGAAAAAAGT') == 1



 
def test_read_write(getSketch):
    fname = str.encode(utils.get_temp_filename('zzz'))
    rng = random.Random(1)
    ctm = getSketch(20, sketchSize)

    kmers = ["".join(rng.choice("ACGT") for _ in range(20))
             for n in range(400)]
    for kmer in kmers:
        ctm.add(kmer)


    ctm.save(fname)

    # on purpose choose parameters that are different from sct
    ctm2 = getSketch.load(fname)
    ctm2.load(fname)
    assert ctm.ksize() == ctm2.ksize()
    for kmer in kmers:
        assert ctm.get(kmer) == ctm2.get(kmer)

        



def test_maxcount_with_bigcount(getSketch):
    # hashtable should not saturate, if use_bigcount is set.
    kh = getSketch(4, 64)

    last_count = None
    for _ in range(0, 10000):
        kh.count('AAAA')
        c = kh.get('AAAA')

        if c == last_count:
            break
        last_count = c

    assert c == 10000, "should be able to count to 1000: %d" % c



# def test_maxcount_with_bigcount_save(getSketch):
#     # hashtable should not saturate, if use_bigcount is set.
#     fname = str.encode(utils.get_temp_filename('zzz'))
#     kh = getSketch(4, 4 ** 4, 4,fname)
#     kh.set_use_bigcount(True)

#     for _ in range(0, 1000):
#         kh.count('AAAA')
#         c = kh.get('AAAA')

#     savepath = utils.get_temp_filename('tempcountingsave.ht')
#     kh.save(savepath)

#     try:
#         kh = Countgraph.load(savepath)
#     except OSError as err:
#         assert 0, "Should not produce an OSError: " + str(err)

#     c = kh.get('AAAA')
#     assert c == 1000, "should be able to count to 1000: %d" % c
#     assert c != MAX_COUNT, c


# def test_bigcount_save():
#     # hashtable should not saturate, if use_bigcount is set.
#     kh = khmer.Countgraph(4, 4 ** 4, 4)
#     kh.set_use_bigcount(True)

#     savepath = utils.get_temp_filename('tempcountingsave.ht')
#     kh.save(savepath)

#     try:
#         kh = Countgraph.load(savepath)
#     except OSError as err:
#         assert 0, "Should not produce an OSError: " + str(err)

#     # set_use_bigcount should still be True after load (i.e. should be saved)

#     assert kh.get('AAAA') == 0

#     for _ in range(0, 1000):
#         kh.count('AAAA')
#         kh.get('AAAA')

#     assert kh.get('AAAA') == 1000


# def test_nobigcount_save():
#     kh = khmer.Countgraph(4, 4 ** 4, 4)
#     # kh.set_use_bigcount(False) <-- this is the default

#     savepath = utils.get_temp_filename('tempcountingsave.ht')
#     kh.save(savepath)

#     try:
#         kh = Countgraph.load(savepath)
#     except OSError as err:
#         assert 0, 'Should not produce an OSError: ' + str(err)

#     # set_use_bigcount should still be False after load (i.e. should be saved)

#     assert kh.get('AAAA') == 0

#     for _ in range(0, 1000):
#         kh.count('AAAA')
#         kh.get('AAAA')

#     assert kh.get('AAAA') == MAX_COUNT


def test_bigcount_abund_dist(getSketch):
    kh =getSketch(18, sketchSize)
    tracking = khmer.Nodegraph(18, 1e2, 4)
    kh.set_use_bigcount(True)

    seqpath = utils.get_test_data('test-abund-read-2.fa')

    kh.consume_seqfile(seqpath)

    dist = kh.abundance_distribution(seqpath, tracking)
    print(kh.get('GGTTGACGGGGCTCAGGG'))

    pdist = [(i, dist[i]) for i in range(len(dist)) if dist[i]]
    assert dist[1002] == 1, pdist



def test_bigcount_overflow(getSketch):
    kh = getSketch(18,1024)
    kh.set_use_bigcount(True)

    for _ in range(0, 65536):
        kh.count('GGTTGACGGGGCTCAGGG')

    assert kh.get('GGTTGACGGGGCTCAGGG') == MAX_BIGCOUNT


def test_get_ksize(getSketch):
    kh = getSketch(22, 16)
    assert kh.ksize() == 22
