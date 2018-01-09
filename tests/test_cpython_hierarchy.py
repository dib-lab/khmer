# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2016, The Regents of the University of California.
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
#     * Neither the name of the University of California nor the names
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
# pylint: disable=C0111,C0103,missing-docstring,no-member,protected-access


import khmer

import pytest
from . import khmer_tst_utils as utils


def test_countgraph_vs_table():
    x = khmer.Counttable(4, 21, 3)
    y = khmer.Countgraph(4, 21, 3)

    assert hasattr(x, 'add')
    assert hasattr(y, 'add')

    assert not hasattr(x, 'consume_and_tag')
    assert hasattr(y, 'consume_and_tag')


def test_nodegraph_vs_table():
    x = khmer.Nodetable(4, 21, 3)
    y = khmer.Nodegraph(4, 21, 3)

    assert hasattr(x, 'add')
    assert hasattr(y, 'add')

    assert not hasattr(x, 'consume_and_tag')
    assert hasattr(y, 'consume_and_tag')


def test_counttable_no_unhash():
    x = khmer.Counttable(4, 21, 3)

    with pytest.raises(ValueError):
        x.reverse_hash(1)


def test_smallcountgraph_vs_table():
    x = khmer.SmallCounttable(4, 21, 3)
    y = khmer.SmallCountgraph(4, 21, 3)

    assert hasattr(x, 'add')
    assert hasattr(y, 'add')

    assert not hasattr(x, 'consume_and_tag')
    assert hasattr(y, 'consume_and_tag')
