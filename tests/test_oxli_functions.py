# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
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
# pylint: disable=missing-docstring,invalid-name,no-member

from . import khmer_tst_utils as utils

# Technically not from 'oxli' but it's what they are
from khmer.khmer_args import (estimate_optimal_with_K_and_M,
                              graphsize_args_report,
                              estimate_optimal_with_K_and_f, optimal_size)


def test_estimate_functions_1():
    res = estimate_optimal_with_K_and_M(99, 1024)
    assert res[0] == 7, res[0]
    assert res[1] == 146, res[1]
    assert res[2] == 1022, res[2]
    assert abs(.008 - res[3]) < .001, res[3]

    res = estimate_optimal_with_K_and_f(99, 0.00701925498897)
    assert res[0] == 7, res[0]
    assert res[1] == 145, res[1]
    assert res[2] == 1015, res[2]
    assert abs(.008 - res[3]) < .002, res[3]

    res = estimate_optimal_with_K_and_M(1024, 2)
    assert res[0] == 1, res[0]
    assert res[1] == 2, res[1]
    assert res[2] == 2, res[2]
    assert res[3] == 1.0, res[3]

    # using a crazy high FP rate just for coverage
    res = estimate_optimal_with_K_and_f(1024, 0.7)
    assert res[0] == 1, res[0]
    assert res[1] == 850, res[1]
    assert res[2] == 850, res[2]
    assert abs(.7 - res[3]) < 0.0022, abs(.7 - res[3])


def test_estimate_functions_namedtup():
    res = estimate_optimal_with_K_and_M(99, 1024)
    assert res.num_htables == 7, res[0]
    assert res.htable_size == 146, res[1]
    assert res.mem_use == 1022, res[2]
    assert abs(.008 - res.fp_rate) < .001, res[3]

    res = estimate_optimal_with_K_and_f(99, 0.00701925498897)
    assert res.num_htables == 7, res[0]
    assert res.htable_size == 145, res[1]
    assert res.mem_use == 1015, res[2]
    assert abs(.008 - res.fp_rate) < .002, res[3]


def test_optimal_size_function():
    res = optimal_size(99, mem_cap=1024)
    assert res.num_htables == 7, res[0]
    assert res.htable_size == 146, res[1]
    assert res.mem_use == 1022, res[2]
    assert abs(.008 - res.fp_rate) < .001, res[3]

    res = optimal_size(99, fp_rate=0.00701925498897)
    assert res.num_htables == 7, res[0]
    assert res.htable_size == 145, res[1]
    assert res.mem_use == 1015, res[2]
    assert abs(.008 - res.fp_rate) < .002, res[3]

    try:
        optimal_size(99, mem_cap=1024, fp_rate=0.00701925498897)
        assert 0, "this should fail"
    except TypeError as err:
        print(str(err))
        assert "num_kmers and either mem_cap or fp_rate" in str(err)

    try:
        optimal_size(99)
        assert 0, "this should fail"
    except TypeError as err:
        print(str(err))
        assert "num_kmers and either mem_cap or fp_rate" in str(err)


def test_output_gen():
    graphsize_args_report(99, 0.00701925498897)
