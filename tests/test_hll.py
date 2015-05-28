#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2014-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# pylint: disable=missing-docstring,protected-access

import math
import string

import khmer

from screed.fasta import fasta_iter

import khmer_tst_utils as utils
from nose.tools import assert_raises


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


def test_hll_empty_fasta():
    filename = utils.get_test_data('test-empty.fa')
    hll = khmer.HLLCounter(ERR_RATE, K)
    with assert_raises(IOError):
        hll.consume_fasta(filename)


def test_hll_consume_fasta():
    # test c++ code to count unique kmers using HyperLogLog

    filename = utils.get_test_data('random-20-a.fa')
    hllcpp = khmer.HLLCounter(ERR_RATE, K)
    hllcpp.consume_fasta(filename)

    assert abs(1 - float(hllcpp.estimate_cardinality()) / N_UNIQUE) < ERR_RATE


def test_hll_consume_fasta_ep():
    # During estimation trigger the _Ep() method,
    # we need all internal counters values to be different than zero for this.

    filename = utils.get_test_data('paired-mixed.fa')
    hll = khmer.HLLCounter(0.36, 32)
    hll.consume_fasta(filename)

    assert all(c != 0 for c in hll.counters)
    assert len(hll) == 236


def test_hll_consume_fasta_estimate_bias():
    # During estimation trigger the estimate_bias method,
    # we need all internal counters values to be different than zero for this,
    # and also the cardinality should be small (if it is large we fall on the
    # default case).

    filename = utils.get_test_data("test-abund-read-3.fa")
    hll = khmer.HLLCounter(0.36, K)
    hll.consume_fasta(filename)

    assert all(c != 0 for c in hll.counters)
    assert len(hll) == 79


def test_hll_len():
    filename = utils.get_test_data('random-20-a.fa')
    hllcpp = khmer.HLLCounter(ERR_RATE, K)
    hllcpp.consume_fasta(filename)

    assert hllcpp.estimate_cardinality() == len(hllcpp)


def test_hll_empty():
    hllcpp = khmer.HLLCounter(ERR_RATE, K)

    assert len(hllcpp) == 0


def test_hll_readonly_alpha():
    hllcpp = khmer.HLLCounter(ERR_RATE, K)
    with assert_raises(AttributeError):
        hllcpp.alpha = 5


def test_hll_cover_calc_alpha():
    hllcpp = khmer.HLLCounter(0.36, K)
    counters = hllcpp.counters
    assert hllcpp.alpha == 0.673
    assert len(counters) == 2 ** 4

    hllcpp = khmer.HLLCounter(0.21, K)
    counters = hllcpp.counters
    assert hllcpp.alpha == 0.697
    assert len(counters) == 2 ** 5

    hllcpp = khmer.HLLCounter(0.16, K)
    counters = hllcpp.counters
    assert hllcpp.alpha == 0.709
    assert len(counters) == 2 ** 6


def test_hll_invalid_base():
    # this test should raise a ValueError,
    # since there are invalid bases in read.

    hllcpp = khmer.HLLCounter(ERR_RATE, 5)
    with assert_raises(ValueError):
        hllcpp.consume_string("ACGTTTCGNAATNNNNN")


def test_hll_invalid_error_rate():
    # test if error_rate is a valid value

    with assert_raises(ValueError):
        hllcpp = khmer.HLLCounter(-0.01, K)


def test_hll_invalid_error_rate_max():
    # test if error_rate is a valid value

    with assert_raises(ValueError):
        hllcpp = khmer.HLLCounter(0.367696, K)


def test_hll_error_rate_max():
    # test if error_rate is a valid value

    hllcpp = khmer.HLLCounter(0.367695, K)
    assert len(hllcpp.counters) == 2 ** 4


def test_hll_invalid_error_rate_min():
    # test if error_rate is a valid value

    with assert_raises(ValueError):
        hllcpp = khmer.HLLCounter(0.0040624, K)


def test_hll_error_rate_min():
    # test if error_rate is a valid value

    hllcpp = khmer.HLLCounter(0.0040625, K)
    assert len(hllcpp.counters) == 2 ** 16


def test_hll_change_error_rate():
    hllcpp = khmer.HLLCounter(0.0040625, K)
    assert hllcpp.error_rate == 0.0040625

    # error rate is discrete, what we test here is if an error rate of 1%
    # rounds to the appropriate value
    hllcpp.error_rate = 0.01
    assert hllcpp.error_rate == 0.008125

    with assert_raises(TypeError):
        del hllcpp.error_rate

    with assert_raises(TypeError):
        hllcpp.error_rate = 5

    with assert_raises(ValueError):
        hllcpp.error_rate = 2.5

    with assert_raises(ValueError):
        hllcpp.error_rate = -10.

    # error rate can only be changed prior to first counting,
    hllcpp.consume_string('AAACCACTTGTGCATGTCAGTGCAGTCAGT')
    with assert_raises(AttributeError):
        hllcpp.error_rate = 0.3


def test_hll_change_ksize():
    hllcpp = khmer.HLLCounter(0.0040625, K)
    assert hllcpp.ksize == K

    hllcpp.ksize = 24
    assert hllcpp.ksize == 24

    hllcpp.ksize = 12L
    assert hllcpp.ksize == 12

    with assert_raises(ValueError):
        hllcpp.ksize = -20

    with assert_raises(TypeError):
        del hllcpp.ksize

    with assert_raises(TypeError):
        hllcpp.ksize = 33.4

    # error rate can only be changed prior to first counting,
    hllcpp.consume_string('AAACCACTTGTGCATGTCAGTGCAGTCAGT')
    with assert_raises(AttributeError):
        hllcpp.ksize = 30


def test_hll_get_counters():
    hll = khmer.HLLCounter(0.36, K)
    counters = hll.counters
    assert len(counters) == 2 ** 4
    assert all(c == 0 for c in counters)
