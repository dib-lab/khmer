#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2014-2015. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=missing-docstring,protected-access

import string

import khmer

from screed.fasta import fasta_iter

import khmer_tst_utils as utils
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


@raises(IOError)
def test_hll_empty_fasta():
    filename = utils.get_test_data('test-empty.fa')
    hll = khmer.HLLCounter(ERR_RATE, K)
    hll.consume_fasta(filename)


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


def test_hll_empty():
    hllcpp = khmer.HLLCounter(ERR_RATE, K)

    assert len(hllcpp) == 0


@raises(AttributeError)
def test_hll_readonly_alpha():
    hllcpp = khmer.HLLCounter(ERR_RATE, K)
    hllcpp.alpha = 5


@raises(AttributeError)
def test_hll_readonly_p():
    hllcpp = khmer.HLLCounter(ERR_RATE, K)
    hllcpp.p = 5


@raises(AttributeError)
def test_hll_readonly_m():
    hllcpp = khmer.HLLCounter(ERR_RATE, K)
    hllcpp.m = 5


def test_hll_cover_calc_alpha():
    hllcpp = khmer.HLLCounter(0.36, K)
    assert hllcpp.alpha == 0.673
    assert hllcpp.p == 4
    assert hllcpp.m == 16

    hllcpp = khmer.HLLCounter(0.21, K)
    assert hllcpp.alpha == 0.697
    assert hllcpp.p == 5
    assert hllcpp.m == 32

    hllcpp = khmer.HLLCounter(0.16, K)
    assert hllcpp.alpha == 0.709
    assert hllcpp.p == 6
    assert hllcpp.m == 64


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

    hllcpp = khmer.HLLCounter(0.367696, K)


def test_hll_error_rate_max():
    # test if error_rate is a valid value

    hllcpp = khmer.HLLCounter(0.367695, K)
    assert hllcpp.p == 4


@raises(ValueError)
def test_hll_invalid_error_rate_min():
    # test if error_rate is a valid value

    hllcpp = khmer.HLLCounter(0.0040624, K)


def test_hll_error_rate_min():
    # test if error_rate is a valid value

    hllcpp = khmer.HLLCounter(0.0040625, K)
    assert hllcpp.p == 16
