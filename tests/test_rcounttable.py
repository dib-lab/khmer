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

from __future__ import print_function

import khmer
import glob

import pytest
from . import khmer_tst_utils as utils


@pytest.mark.parametrize('sketchtype', [
    (khmer.RNodetable),
    (khmer.RCounttable),
    (khmer.RSmallCounttable),
])
def test_basic(sketchtype):
    """
    Basic sanity check

    Loading a sequence with 3 k-mers should give the same result as loading the
    k-mers separately.
    """
    sketch1 = sketchtype(25, 1e4, 4)
    allseq = utils.get_test_data('rollhash-allseq.fa')
    sketch1.consume(allseq)

    sketch2 = sketchtype(25, 1e4, 4)
    for f in glob.glob('rollhash-seq?.fa'):
        seqfile = utils.get_test_data(f)
        sketch2.consume(seqfile)

    assert sketch1.get('AAGATTTCCCGCCCCGATCGATTCG') > 0
    assert sketch1.get('AGATTTCCCGCCCCGATCGATTCGT') > 0
    assert sketch1.get('GATTTCCCGCCCCGATCGATTCGTT') > 0

    assert sketch2.get('AAGATTTCCCGCCCCGATCGATTCG') > 0
    assert sketch2.get('AGATTTCCCGCCCCGATCGATTCGT') > 0
    assert sketch2.get('GATTTCCCGCCCCGATCGATTCGTT') > 0
