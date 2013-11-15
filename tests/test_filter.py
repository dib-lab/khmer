#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import khmer
from screed.fasta import fasta_iter
from nose.plugins.attrib import attr

import khmer_tst_utils as utils


def teardown():
    utils.cleanup()


def load_fa_seq_names(filename):
    fp = open(filename)
    records = list(fasta_iter(fp))
    names = [r['name'] for r in records]
    return names


class Test_Filter(object):

    @attr('highmem')
    def test_abund(self):
        ht = khmer.new_hashtable(10, 4 ** 10)

        filename = utils.get_test_data('test-abund-read.fa')
        outname = utils.get_temp_filename('test_abund.out')

        ht.consume_fasta(filename)
        ht.output_fasta_kmer_pos_freq(filename, outname)

        fd = open(outname, "r")

        output = fd.readlines()
        assert len(output) == 1

        output = output[0]
        output = output.strip().split()

        assert ['1'] * (114 - 10 + 1) == output

        fd.close()


@attr('highmem')
def test_filter_sodd():
    K = 32
    HASHTABLE_SIZE = int(8e7)
    N_HT = 4
    MAX_SODD = 3

    ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)
    filename = utils.get_test_data('../../data/high-sodd.fa')

    ht.consume_fasta(filename)

    seq = "CGTTAGTTGCGGTGCCGACCGGCAAACTTGGTTTTGCCAAAAATTTTTACAGTTAGAAATTATTC" \
          "ACAAAGTTGCACCGGAATTCGGTTACAAACGTCATTCTAACTAAT"
    trim_seq, trim_at = ht.trim_on_sodd(seq, MAX_SODD)
    assert trim_seq == "CGTTAGTTGCGGTGCCGACCGGCAAACTTGGT"

    seq = "ACAAAATTCCACATATAGTCATAATTGTGGGCAATTTTCGTCCCAAATTAGTTAGAATGACGTTT" \
          "GTAACCGAATTCCGGTGCAACTTTGTGAATAATTTCTAACTGTAAAAAT"
    trim_seq, trim_at = ht.trim_on_sodd(seq, MAX_SODD)
    assert trim_seq == "ACAAAATTCCACATATAGTCATAATTGTGGGCAATT"

    seq = "GCACGCAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTG"
    trim_seq, trim_at = ht.trim_on_sodd(seq, MAX_SODD)
    assert trim_seq == seq
