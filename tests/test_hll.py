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
from nose.tools import raises


TT = string.maketrans('ACGT', 'TGCA')


def teardown():
    utils.cleanup()


def test_hll_add_python():
    # test python code to count unique kmers using HyperLogLog.
    # use the lower level add() method, which accepts anything,
    # and compare to an exact count using collections.Counter

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


def test_hll_consume_string():
    # test c++ code to count unique kmers using HyperLogLog,
    # using screed to feed each read to the counter.

    filename = utils.get_test_data('random-20-a.fa')

    K = 20  # size of kmer
    ERROR_RATE = 0.01

    hllcpp = khmer.new_hll_counter(ERROR_RATE, K)

    for n, record in enumerate(fasta_iter(open(filename))):
        hllcpp.consume_string(record['sequence'])

    assert abs(1 - float(hllcpp.estimate_cardinality()) / 3960) < 0.01


def test_hll_consume_fasta():
    # test c++ code to count unique kmers using HyperLogLog

    filename = utils.get_test_data('random-20-a.fa')

    K = 20  # size of kmer
    ERROR_RATE = 0.01

    hllcpp = khmer.new_hll_counter(ERROR_RATE, K)

    hllcpp.consume_fasta(filename)

    assert abs(1 - float(hllcpp.estimate_cardinality()) / 3960) < 0.01


@raises(ValueError)
def test_hll_invalid_base():
    # this test should raise a ValueError,
    # since there are invalid bases in read.

    K = 5  # size of kmer
    ERROR_RATE = 0.01

    hllcpp = khmer.new_hll_counter(ERROR_RATE, K)

    hllcpp.consume_string("ACGTTTCGNAATNNNNN")


@raises(ValueError)
def test_hll_invalid_error_rate():
    # test if error_rate is a valid value

    K = 20  # size of kmer
    ERROR_RATE = -0.01

    hllcpp = khmer.new_hll_counter(ERROR_RATE, K)


@raises(ValueError)
def test_hll_invalid_error_rate_max():
    # test if error_rate is a valid value

    K = 20  # size of kmer
    ERROR_RATE = 0.50

    hllcpp = khmer.new_hll_counter(ERROR_RATE, K)


@raises(ValueError)
def test_hll_invalid_error_rate_min():
    # test if error_rate is a valid value

    K = 20  # size of kmer
    ERROR_RATE = 0.000001

    hllcpp = khmer.HLLCounter(ERROR_RATE, K)
