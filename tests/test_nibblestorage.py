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
# Contact: khmer-project@idyll.org


import random

from khmer import SmallCounttable

from . import khmer_tst_utils as utils


# This mostly tests the underlying NibbleStorage class

def test_single_add():
    sct = SmallCounttable(4, 1e6, 4)
    sct.add("AAAA")
    assert sct.get("AAAA") == 1


def test_split_byte_murmur():
    # check the byte is correctly split when using murmur hash
    sct = SmallCounttable(4, 4, 1)

    # these kmers were carefully chosen to have hash values that produce
    # consecutive indices in the count table.
    a = "AAAC"
    b = "AAAG"

    assert sct.get_kmer_hashes(a) == [11898086063751343884]
    assert sct.get_kmer_hashes(b) == [10548630838975263317]

    sct.add(a)

    assert sct.get(a) == 1
    assert sct.get(b) == 0


def test_overflow():
    # check that we do not overflow into other part of the byte
    sct = SmallCounttable(4, 1e6, 4)
    a = "AAAA"
    b = "AAAT"

    # try to overflow our 4bit counter
    for n in range(17):
        sct.add(a)

    assert sct.get(a) == 15
    assert sct.get(b) == 0

    # repeat with the other kmer that hashes to the other half of the byte
    sct = SmallCounttable(4, 1e6, 4)
    a = "AAAA"
    b = "AAAT"

    # try to overflow our 4bit counter
    for n in range(17):
        sct.add(b)

    assert sct.get(b) == 15
    assert sct.get(a) == 0


def test_random_kmers():
    # check for out-of-bounds errors and similar with random kmers
    rng = random.Random(1)

    sct = SmallCounttable(20, 1e2, 4)

    kmers = ["".join(rng.choice("ACGT") for _ in range(20))
             for n in range(400)]
    for kmer in kmers:
        sct.add(kmer)

    for kmer in kmers:
        sct.get(kmer)


def test_read_write():
    rng = random.Random(1)

    sct = SmallCounttable(20, 1e2, 4)

    kmers = ["".join(rng.choice("ACGT") for _ in range(20))
             for n in range(400)]
    for kmer in kmers:
        sct.add(kmer)

    fname = utils.get_temp_filename('zzz')

    sct.save(fname)

    # on purpose choose parameters that are different from sct
    sct2 = SmallCounttable.load(fname)
    assert sct.ksize() == sct2.ksize()
    for kmer in kmers:
        assert sct.get(kmer) == sct2.get(kmer)
