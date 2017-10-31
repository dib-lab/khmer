# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2013-2015, Michigan State University.
# Copyright (C) 2015-2016, The Regents of the University of California.
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
# pylint: disable=missing-docstring,invalid-name


from khmer import Countgraph, SmallCountgraph, Nodegraph
from khmer import (Nodetable, Counttable, CyclicCounttable, SmallCounttable,
                   QFCounttable)
from khmer._oxli.utils import get_n_primes_near_x

import math
import pytest

params_1m = (1000003, 2)
PRIMES_1m = [100003, 1000007]
QF_SIZE = 2**math.ceil(math.log(PRIMES_1m[0], 2))


def tablewrapper(tabletype):

    def build(k, *args):
        try:
            starting_size, n_tables = args
        except:
            starting_size, n_tables = params_1m

        if tabletype is QFCounttable:
            qf_size = 2**math.ceil(math.log(starting_size, 2))
            return tabletype(k, qf_size)
        else:
            return tabletype(k, starting_size, n_tables)

    return build


@pytest.fixture(params=[Countgraph, Counttable, CyclicCounttable,
                        SmallCountgraph, SmallCounttable, Nodegraph,
                        Nodetable])
def Tabletype(request):
    return tablewrapper(request.param)


# all the table types!
@pytest.fixture(params=[Countgraph, Counttable, SmallCountgraph,
                        SmallCounttable, Nodegraph, Nodetable,
                        QFCounttable])
def AnyTabletype(request):
    return tablewrapper(request.param)


# all the counting types!
@pytest.fixture(params=[Countgraph, Counttable, CyclicCounttable,
                        SmallCountgraph, SmallCounttable])
def Countingtype(request):
    return tablewrapper(request.param)


# all the graph types!
@pytest.fixture(params=[Countgraph, Nodegraph])
def Graphtype(request):
    return tablewrapper(request.param)
