from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#

# pylint: disable=C0111,C0103,E1103,W0612

from . import khmer_tst_utils as utils
import khmer
from oxli import functions


def test_estimate_functions_1():
    res = functions.estimate_optimal_with_K_and_M(99, 1024)
    assert res[0] == 7, res[0]
    assert res[1] == 146, res[1]
    assert res[2] == 1022, res[2]
    assert abs(.008 - res[3]) < .001, res[3]

    res = functions.estimate_optimal_with_K_and_f(99, 0.00701925498897)
    assert res[0] == 7, res[0]
    assert res[1] == 145, res[1]
    assert res[2] == 1015, res[2]
    assert abs(.008 - res[3]) < .002, res[3]

    res = functions.estimate_optimal_with_K_and_M(1024, 2)
    assert res[0] == 1, res[0]
    assert res[1] == 2, res[1]
    assert res[2] == 2, res[2]
    assert res[3] == 1.0, res[3]

    # using a crazy high FP rate just for coverage
    res = functions.estimate_optimal_with_K_and_f(1024, 0.7)
    assert res[0] == 1, res[0]
    assert res[1] == 850, res[1]
    assert res[2] == 850, res[2]
    assert abs(.7 - res[3]) < 0.0022, abs(.7 - res[3])


def test_estimate_functions_namedtup():
    res = functions.estimate_optimal_with_K_and_M(99, 1024)
    assert res.num_htables == 7, res[0]
    assert res.htable_size == 146, res[1]
    assert res.mem_use == 1022, res[2]
    assert abs(.008 - res.fp_rate) < .001, res[3]

    res = functions.estimate_optimal_with_K_and_f(99, 0.00701925498897)
    assert res.num_htables == 7, res[0]
    assert res.htable_size == 145, res[1]
    assert res.mem_use == 1015, res[2]
    assert abs(.008 - res.fp_rate) < .002, res[3]


def test_optimal_size_function():
    res = functions.optimal_size(99, mem_cap=1024)
    assert res.num_htables == 7, res[0]
    assert res.htable_size == 146, res[1]
    assert res.mem_use == 1022, res[2]
    assert abs(.008 - res.fp_rate) < .001, res[3]

    res = functions.optimal_size(99, fp_rate=0.00701925498897)
    assert res.num_htables == 7, res[0]
    assert res.htable_size == 145, res[1]
    assert res.mem_use == 1015, res[2]
    assert abs(.008 - res.fp_rate) < .002, res[3]

    try:
        functions.optimal_size(99, mem_cap=1024, fp_rate=0.00701925498897)
        assert 0, "this should fail"
    except TypeError as err:
        print(str(err))
        assert "num_kmers and either mem_cap or fp_rate" in str(err)

    try:
        functions.optimal_size(99)
        assert 0, "this should fail"
    except TypeError as err:
        print(str(err))
        assert "num_kmers and either mem_cap or fp_rate" in str(err)


def test_output_gen():
    res = functions.optimal_args_output_gen(99, 0.00701925498897)
