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
# pylint: disable=missing-docstring,protected-access

from __future__ import print_function
from __future__ import absolute_import, unicode_literals

import khmer
from . import khmer_tst_utils as utils
import screed

# add:
# * get default params from Python
# * keyword args for minhash constructor
# * trap error from handing protein/non-DNA to a DNA MH
# * fail on untagged/unloaded countgraph
# * nan on empty minhash

def test_default_params():
    # verify that MHs have these default parameters.
    mh = khmer.MinHash(100, 32)
    mh2 = khmer.MinHash(100, 32, 9999999967, False)  # should not be changed!

    mh.add_hash(9999999968)
    mh2.add_hash(9999999968)
    assert mh.get_mins() == mh2.get_mins()
    assert mh.get_mins()[0] == 1

def test_basic_dna():
    # verify that MHs of size 1 stay size 1, & act properly as bottom sketches.
    mh = khmer.MinHash(1, 4)
    mh.add_sequence('ATGC')
    a = mh.get_mins()

    mh.add_sequence('GCAT')             # this will not get added; hash > ATGC
    b = mh.get_mins()

    print(a, b)
    assert a == b
    assert len(b) == 1


def test_protein():
    # verify that we can hash protein/aa sequences
    mh = khmer.MinHash(1, 4, 9999999967, True)
    mh.add_protein('AGYYG')

    assert len(mh.get_mins()) == 1


def test_build_nbhd():
    # load in a sequence,
    # calculate neighborhood minhashes,
    # combine nbhd into combined minhashes,
    # verify that there is only one.
    ct = khmer.Countgraph(32, 1e7, 4)
    inpath = utils.get_test_data('2kb-random.fa')

    for record in screed.open(inpath):
        ct.consume_and_tag(record.sequence)

    nbhd_mh = ct.build_neighborhood_minhashes(20)
    combined = nbhd_mh.build_combined_minhashes(500)

    assert len(combined) == 1, combined


def test_build_save_load_nbhd():
    # load in a sequence,
    # calculate neighborhood minhashes,
    # save, load
    # combine nbhd into combined minhashes,
    # verify that there is only one.
    ct = khmer.Countgraph(32, 1e7, 4)
    inpath = utils.get_test_data('2kb-random.fa')
    savepath = utils.get_temp_filename('save.mhi')

    for record in screed.open(inpath):
        ct.consume_and_tag(record.sequence)

    nbhd_mh = ct.build_neighborhood_minhashes(20)
    nbhd_mh.save(savepath)
    del nbhd_mh

    nbhd_mh = khmer.load_neighborhood_minhash(savepath)
    combined = nbhd_mh.build_combined_minhashes(500)

    assert len(combined) == 1, combined

def test_build_nbhd_2():
    # load in two disconnected sequences,
    # calculate neighborhood minhashes,
    # combine nbhd into combined minhashes,
    # verify that there are two.
    ct = khmer.Countgraph(32, 1e7, 4)
    inpath = utils.get_test_data('2kb-random.fa')
    inpath2 = utils.get_test_data('2kb-random-b.fa')

    for record in screed.open(inpath):
        ct.consume_and_tag(record.sequence)
    for record in screed.open(inpath2):
        ct.consume_and_tag(record.sequence)

    nbhd_mh = ct.build_neighborhood_minhashes(20)
    combined = nbhd_mh.build_combined_minhashes(500)

    assert len(combined) == 2, combined


def test_build_save_load_nbhd_2():
    # load in two disconnected sequences,
    # calculate neighborhood minhashes,
    # save, load
    # combine nbhd into combined minhashes,
    # verify that there are two.
    ct = khmer.Countgraph(32, 1e7, 4)
    inpath = utils.get_test_data('2kb-random.fa')
    inpath2 = utils.get_test_data('2kb-random-b.fa')
    savepath = utils.get_temp_filename('save.mhi')

    for record in screed.open(inpath):
        ct.consume_and_tag(record.sequence)
    for record in screed.open(inpath2):
        ct.consume_and_tag(record.sequence)

    nbhd_mh = ct.build_neighborhood_minhashes(20)
    nbhd_mh.save(savepath)
    del nbhd_mh

    nbhd_mh = khmer.load_neighborhood_minhash(savepath)
    combined = nbhd_mh.build_combined_minhashes(500)

    assert len(combined) == 2, combined

def test_build_combined_minhashes_2():
    # load in two (disconnected) sequences,
    # calculate neighborhood minhashes,
    # combine them into combined minhashes,
    # verify that the signatures of the combined minhashes match the
    #    minhash signatures of the two input sequences

    ct = khmer.Countgraph(32, 1e7, 4)
    inpath = utils.get_test_data('2kb-random.fa')
    inpath2 = utils.get_test_data('2kb-random-b.fa')

    mh1 = khmer.MinHash(1000, 32)
    for record in screed.open(inpath):
        ct.consume_and_tag(record.sequence)
        mh1.add_sequence(record.sequence)
        
    mh2 = khmer.MinHash(1000, 32)
    for record in screed.open(inpath2):
        ct.consume_and_tag(record.sequence)
        mh2.add_sequence(record.sequence)

    nbhd_mh = ct.build_neighborhood_minhashes(20)
    combined = nbhd_mh.build_combined_minhashes(1000)

    assert len(combined) == 2

    # do we find the first signature in combined?
    found = False
    for c in combined:
        print(mh1.compare(c.get_minhash()))
        if mh1.compare(c.get_minhash()) > 0.4:   # why not ~1?? CTB
            combined.remove(c)
            found = True
    assert found, "didn't find 2kb-random signature in combined"

    # now that we've removed that match, do we find another? shouldn't.
    found = False
    for c in combined:
        print(mh1.compare(c.get_minhash()))
        if mh1.compare(c.get_minhash()) > 0.4:   # why not ~1?? CTB
            found = True
    assert not found, "found the 2kb-random signature that I removed!"

    # make sure we find the second signature, too.
    found = False
    for c in combined:
        print(mh2.compare(c.get_minhash()))
        if mh2.compare(c.get_minhash()) > 0.4:   # why not ~1?? CTB
            found = True
    assert found, "didn't find 2kb-random-b signature in combined"
