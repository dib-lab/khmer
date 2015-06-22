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
    res = functions.estimate_optimal_with_N_and_M(99, 1024)
    assert res[0] == 7, res[0]
    assert res[1] == 146, res[1]
    assert res[2] == 1022, res[2]
    assert abs(.008 - res[3]) < .001, res[3]

    res = functions.estimate_optimal_with_N_and_f(99, 0.00701925498897)
    assert res[0] == 7, res[0]
    assert res[1] == 145, res[1]
    assert res[2] == 1015, res[2]
    assert abs(.008 - res[3]) < .002, res[3]

    res = functions.estimate_optimal_with_N_and_M(1024, 2)
    assert res[0] == 1, res[0]
    assert res[1] == 2, res[1]
    assert res[2] == 2, res[2]
    assert res[3] == 1.0, res[3]

    # using a crazy high FP rate just for coverage
    res = functions.estimate_optimal_with_N_and_f(1024, 0.7)
    assert res[0] == 1, res[0]
    assert res[1] == 850, res[1]
    assert res[2] == 850, res[2]
    assert abs(.7 - res[3]) < 0.0022, abs(.7 - res[3])
