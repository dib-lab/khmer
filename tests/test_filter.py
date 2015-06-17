from __future__ import print_function
from __future__ import absolute_import
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
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
        ht = khmer.new_hashtable(10, 4 ** 10)

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
        except IOError as err:
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
