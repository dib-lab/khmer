#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=missing-docstring,protected-access

import string

import khmer
from khmer import ReadParser

from screed.fasta import fasta_iter
import screed

import khmer_tst_utils as utils
from nose.plugins.attrib import attr


TT = string.maketrans('ACGT', 'TGCA')


def teardown():
    utils.cleanup()


def test_hll_python_1():
    # test python code to count unique kmers using HyperLogLog
    filename = utils.get_test_data('random-20-a.fa')

    K = 20  # size of kmer
    ERROR_RATE = 0.01  # bounded error

    hllcpp = khmer.new_hll_counter(ERROR_RATE, K)
    counter = set()

    for n, record in enumerate(fasta_iter(open(filename))):
        sequence = record['sequence']
        seq_len = len(sequence)
        for n in range(0, seq_len + 1 - K):
            kmer = sequence[n:n + K]
            rc = kmer[::-1].translate(TT)

            hllcpp.add(kmer)

            if rc in counter:
                kmer = rc
            counter.update([kmer])

    n_unique = len(counter)

    assert n_unique == 3960
    assert abs(1 - float(hllcpp.estimate_cardinality()) / n_unique) < 0.01


def test_hll_c_1():
    # test c++ code to count unique kmers using HyperLogLog

    filename = utils.get_test_data('random-20-a.fa')

    K = 20  # size of kmer
    ERROR_RATE = 0.01

    hllcpp = khmer.new_hll_counter(ERROR_RATE, K)

    for n, record in enumerate(fasta_iter(open(filename))):
        hllcpp.consume_string(record['sequence'])

    assert abs(1 - float(hllcpp.estimate_cardinality()) / 3960) < 0.01
