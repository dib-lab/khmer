#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2014-2015. It is licensed under
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
from nose.tools import raises


TT = string.maketrans('ACGT', 'TGCA')
K = 20  # size of kmer
ERR_RATE = 0.01
N_UNIQUE = 3960


def teardown():
    utils.cleanup()


def test_hll_add_python():
    # test python code to count unique kmers using HyperLogLog.
    # use the lower level add() method, which accepts anything,
    # and compare to an exact count using collections.Counter

    filename = utils.get_test_data('random-20-a.fa')
    hllcpp = khmer.HLLCounter(ERR_RATE, K)
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

    assert n_unique == N_UNIQUE
    assert abs(1 - float(hllcpp.estimate_cardinality()) / N_UNIQUE) < ERR_RATE


def test_hll_consume_string():
    # test c++ code to count unique kmers using HyperLogLog,
    # using screed to feed each read to the counter.

    filename = utils.get_test_data('random-20-a.fa')
    hllcpp = khmer.HLLCounter(ERR_RATE, K)
    for n, record in enumerate(fasta_iter(open(filename))):
        hllcpp.consume_string(record['sequence'])

    assert abs(1 - float(hllcpp.estimate_cardinality()) / N_UNIQUE) < ERR_RATE


def test_hll_consume_fasta():
    # test c++ code to count unique kmers using HyperLogLog

    filename = utils.get_test_data('random-20-a.fa')
    hllcpp = khmer.HLLCounter(ERR_RATE, K)
    hllcpp.consume_fasta(filename)

    assert abs(1 - float(hllcpp.estimate_cardinality()) / N_UNIQUE) < ERR_RATE


def test_hll_len():
    filename = utils.get_test_data('random-20-a.fa')
    hllcpp = khmer.HLLCounter(ERR_RATE, K)
    hllcpp.consume_fasta(filename)

    assert hllcpp.estimate_cardinality() == len(hllcpp)


def test_hll_merge():
    hll_total = khmer.HLLCounter(ERR_RATE, K)
    hll_merged = khmer.HLLCounter(ERR_RATE, K)

    filename = utils.get_test_data("test-abund-read-2.fa")
    hll_partial_1 = khmer.HLLCounter(ERR_RATE, K)
    hll_partial_1.consume_fasta(filename)
    hll_total.consume_fasta(filename)

    filename = utils.get_test_data("test-abund-read-3.fa")
    hll_partial_2 = khmer.HLLCounter(ERR_RATE, K)
    hll_partial_2.consume_fasta(filename)
    hll_total.consume_fasta(filename)

    hll_merged.merge(hll_partial_1)
    hll_merged.merge(hll_partial_2)

    assert len(hll_total) == len(hll_merged)


@raises(ValueError)
def test_hll_invalid_base():
    # this test should raise a ValueError,
    # since there are invalid bases in read.

    hllcpp = khmer.HLLCounter(ERR_RATE, 5)
    hllcpp.consume_string("ACGTTTCGNAATNNNNN")


@raises(ValueError)
def test_hll_invalid_error_rate():
    # test if error_rate is a valid value

    hllcpp = khmer.HLLCounter(-0.01, K)


@raises(ValueError)
def test_hll_invalid_error_rate_max():
    # test if error_rate is a valid value

    hllcpp = khmer.HLLCounter(0.5, K)


@raises(ValueError)
def test_hll_invalid_error_rate_min():
    # test if error_rate is a valid value

    hllcpp = khmer.HLLCounter(0.000001, K)
