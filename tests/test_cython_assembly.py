# -*- coding: UTF-8 -*-
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2010-2015, Michigan State University.
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
# pylint: disable=missing-docstring,protected-access,no-member,invalid-name

from __future__ import print_function
from __future__ import absolute_import

import itertools
import random

import khmer
from khmer.khmer_args import estimate_optimal_with_K_and_f as optimal_fp
from khmer import ReadParser
from khmer import reverse_complement as revcomp
from . import khmer_tst_utils as utils
from khmer._oxli.assembly import LinearAssembler

import pytest
import screed

from .graph_features import *
from .graph_features import K

def teardown():
    utils.cleanup()

@pytest.mark.parametrize("assembler", [LinearAssembler])
class TestNonBranching:

    def test_all_start_positions(self, linear_structure, assembler):
        # assemble entire contig, starting from wherever
        graph, contig = linear_structure
        asm = assembler(graph)

        for start in range(0, len(contig), 150):
            path = asm.assemble(contig[start:start + K])
            assert utils._equals_rc(path, contig), start

    def test_all_left_to_beginning(self, linear_structure, assembler):
        # assemble directed left
        graph, contig = linear_structure
        asm = assembler(graph)

        for start in range(0, len(contig), 150):
            path = asm.assemble_left(contig[start:start + K])
            print(path, ', ', contig[:start])
            assert utils._equals_rc(path, contig[:start+K]), start

    def test_all_right_to_end(self, linear_structure, assembler):
        # assemble directed right
        graph, contig = linear_structure
        asm = assembler(graph)

        for start in range(0, len(contig), 150):
            path = asm.assemble_right(contig[start:start + K])
            print(path, ', ', contig[:start])
            assert utils._equals_rc(path, contig[start:]), start


