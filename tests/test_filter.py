# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2010-2015, Michigan State University.
# Copyright (C) 2015, The Regents of the University of California.
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
from __future__ import print_function
from __future__ import absolute_import
import khmer
from screed.fasta import fasta_iter
from nose.plugins.attrib import attr

from . import khmer_tst_utils as utils


def teardown():
    utils.cleanup()


def load_fa_seq_names(filename):
    fp = open(filename)
    records = list(fasta_iter(fp))
    names = [r['name'] for r in records]
    return names


class Test_Filter(object):

    def test_abund(self):
        ht = khmer.Countgraph(10, 4 ** 10, 1)

        filename = utils.get_test_data('test-abund-read.fa')
        outname = utils.get_temp_filename('test_abund.out')

        ht.consume_fasta(filename)
        try:
            ht.consume_fasta()
            assert 0, "should fail"
        except TypeError as err:
            print(str(err))
        try:
            ht.consume_fasta("nonexistent")
            assert 0, "should fail"
        except OSError as err:
            print(str(err))
        ht.output_fasta_kmer_pos_freq(filename, outname)
        try:
            ht.output_fasta_kmer_pos_freq()
            assert 0, "should fail"
        except TypeError as err:
            print(str(err))

        fd = open(outname, "r")

        output = fd.readlines()
        assert len(output) == 1

        output = output[0]
        output = output.strip().split()

        assert ['1'] * (114 - 10 + 1) == output

        fd.close()
