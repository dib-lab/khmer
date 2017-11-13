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
            assert utils._equals_rc(path, contig[:start + K]), start

    def test_all_right_to_end(self, linear_structure, assembler):
        # assemble directed right
        graph, contig = linear_structure
        asm = assembler(graph)

        for start in range(0, len(contig), 150):
            path = asm.assemble_right(contig[start:start + K])
            print(path, ', ', contig[:start])
            assert utils._equals_rc(path, contig[start:]), start

    def test_circular(self, circular_linear_structure, assembler):

        graph, contig = circular_linear_structure
        asm = assembler(graph)

        path = asm.assemble_right(contig[:K])
        print(path, ',', contig)
        assert utils._equals_rc(path, contig[:len(path)])

    def test_hash_as_seed(self, linear_structure, assembler):
        graph, contig = linear_structure
        asm = assembler(graph)

        left = graph.hash(contig[:K])
        assert utils._equals_rc(asm.assemble(left), contig)


class TestLinearAssembler_RightBranching:

    def test_branch_point(self, right_tip_structure):
        graph, contig, L, HDN, R, tip = right_tip_structure

        assert graph.kmer_degree(HDN) == 3

    def test_beginning_to_branch(self, right_tip_structure):
        # assemble from beginning of contig, up until branch point
        graph, contig, L, HDN, R, tip = right_tip_structure
        asm = khmer.LinearAssembler(graph)
        path = asm.assemble(contig[0:K])

        assert len(path) == HDN.pos + K
        assert utils._equals_rc(path, contig[:len(path)])

    def test_assemble_takes_hash(self, right_tip_structure):
        # assemble from beginning of contig, up until branch point
        graph, contig, L, HDN, R, tip = right_tip_structure
        asm = khmer.LinearAssembler(graph)
        path = asm.assemble(graph.hash(contig[0:K]))

        assert len(path) == HDN.pos + K
        assert utils._equals_rc(path, contig[:len(path)])

    def test_beginning_to_branch_revcomp(self, right_tip_structure):
        # assemble from beginning of contig, up until branch point
        # starting from rev comp
        graph, contig, L, HDN, R, tip = right_tip_structure
        asm = khmer.LinearAssembler(graph)
        path = asm.assemble(revcomp(contig[0:K]))

        assert len(path) == HDN.pos + K
        assert utils._equals_rc(path, contig[:len(path)])

    def test_left_of_branch_to_beginning(self, right_tip_structure):
        # start from HDN (left of branch)
        graph, contig, L, HDN, R, tip = right_tip_structure
        asm = khmer.LinearAssembler(graph)
        path = asm.assemble(L)

        assert len(path) == HDN.pos + K
        assert utils._equals_rc(path, contig[:len(path)])

    def test_left_of_branch_to_beginning_revcomp(self, right_tip_structure):
        # start from revcomp of HDN (left of branch)
        graph, contig, L, HDN, R, tip = right_tip_structure
        asm = khmer.LinearAssembler(graph)
        path = asm.assemble(revcomp(L))

        assert len(path) == HDN.pos + K
        assert utils._equals_rc(path, contig[:len(path)])

    def test_right_of_branch_outwards_to_ends(self, right_tip_structure):
        # assemble from right of branch point (at R)
        # Should get the *entire* original contig, as the assembler
        # will move left relative to the branch, and not consider it
        # as a high degree node
        graph, contig, L, HDN, R, tip = right_tip_structure
        asm = khmer.LinearAssembler(graph)
        path = asm.assemble(R)

        assert len(path) == len(contig)
        assert utils._equals_rc(path, contig)

    def test_end_to_beginning(self, right_tip_structure):
        # should have exact same behavior as right_of_branch_outwards
        graph, contig, L, HDN, R, tip = right_tip_structure
        asm = khmer.LinearAssembler(graph)
        path = asm.assemble(contig[-K:])

        assert len(path) == len(contig)
        assert utils._equals_rc(path, contig)


class TestLinearAssembler_LeftBranching:

    def test_branch_point(self, left_tip_structure):
        graph, contig, L, HDN, R, tip = left_tip_structure

        assert graph.kmer_degree(HDN) == 3

    def test_end_to_branch(self, left_tip_structure):
        # assemble from end until branch point
        # should include HDN
        graph, contig, L, HDN, R, tip = left_tip_structure
        asm = khmer.LinearAssembler(graph)
        path = asm.assemble(contig[-K:])

        assert len(path) == len(contig) - HDN.pos
        assert utils._equals_rc(path, contig[HDN.pos:])

    def test_branch_to_end(self, left_tip_structure):
        # assemble from branch point until end
        graph, contig, L, HDN, R, tip = left_tip_structure
        asm = khmer.LinearAssembler(graph)
        path = asm.assemble(HDN)

        assert len(path) == len(contig) - HDN.pos
        assert utils._equals_rc(path, contig[HDN.pos:])

    def test_from_branch_to_ends_with_stopbf(self, left_tip_structure):
        # block the tip with the stop_filter. should return a full length
        # contig.
        graph, contig, L, HDN, R, tip = left_tip_structure

        stop_filter = khmer.Nodegraph(K, 1e5, 4)
        stop_filter.count(tip)

        asm = khmer.LinearAssembler(graph, stop_filter=stop_filter)

        path = asm.assemble(HDN)

        assert len(path) == len(contig)
        assert utils._equals_rc(path, contig)

    def test_from_branch_to_ends_with_stopbf_revcomp(self, left_tip_structure):
        # block the tip with the stop_filter. should return a full length
        # contig.
        graph, contig, L, HDN, R, tip = left_tip_structure

        stop_filter = khmer.Nodegraph(K, 1e5, 4)
        stop_filter.count(tip)
        asm = khmer.LinearAssembler(graph, stop_filter=stop_filter)

        path = asm.assemble(revcomp(HDN))

        assert len(path) == len(contig)
        assert utils._equals_rc(path, contig)

    def test_end_thru_tip_with_stopbf(self, left_tip_structure):
        # assemble up to branch point, and include introduced branch b/c
        # of stop bf
        graph, contig, L, HDN, R, tip = left_tip_structure

        stop_filter = khmer.Nodegraph(K, 1e5, 4)
        stop_filter.count(L)          # ...and block original path
        asm = khmer.LinearAssembler(graph, stop_filter=stop_filter)

        path = asm.assemble(contig[-K:])
        assert len(path) == len(contig) - HDN.pos + 1

        # should be the tip k-kmer, plus the last base of the HDN thru
        # the end of the contig
        assert utils._equals_rc(path, tip + contig[HDN.pos + K - 1:])

    def test_single_node_flanked_by_hdns(self, left_tip_structure):
        # assemble single node flanked by high-degree nodes
        # we'll copy the main nodegraph before mutating it
        graph, contig, L, HDN, R, tip = left_tip_structure
        asm = khmer.LinearAssembler(graph)

        graph.consume(mutate_position(contig, HDN.pos + K))

        path = asm.assemble(HDN)

        assert len(path) == K
        assert utils._equals_rc(path, HDN)


class TestLabeledAssembler:

    def test_hash_as_seed(self, linear_structure):
        graph, contig = linear_structure
        lh = khmer.GraphLabels(graph)
        asm = khmer.SimpleLabeledAssembler(lh)

        left = graph.hash(contig[:K])
        assert utils._equals_rc(asm.assemble(left).pop(), contig)

    def test_beginning_to_end_across_tip(self, right_tip_structure):
        # assemble entire contig, ignoring branch point b/c of labels
        graph, contig, L, HDN, R, tip = right_tip_structure
        lh = khmer.GraphLabels(graph)
        asm = khmer.SimpleLabeledAssembler(lh)
        hdn = graph.find_high_degree_nodes(contig)
        # L, HDN, and R will be labeled with 1
        lh.label_across_high_degree_nodes(contig, hdn, 1)

        path = asm.assemble(contig[:K])

        assert len(path) == 1, "there should only be one path"
        path = path[0]  # @CTB

        assert len(path) == len(contig)
        assert utils._equals_rc(path, contig)

    def test_assemble_right_double_fork(self, right_double_fork_structure):
        # assemble two contigs from a double forked structure
        graph, contig, L, HDN, R, branch = right_double_fork_structure
        lh = khmer.GraphLabels(graph)
        asm = khmer.SimpleLabeledAssembler(lh)

        hdn = graph.find_high_degree_nodes(contig)
        hdn += graph.find_high_degree_nodes(branch)
        print(list(hdn))
        lh.label_across_high_degree_nodes(contig, hdn, 1)
        lh.label_across_high_degree_nodes(branch, hdn, 2)
        print(lh.get_tag_labels(list(hdn)[0]))

        paths = asm.assemble(contig[:K])
        print('Path lengths', [len(x) for x in paths])

        assert len(paths) == 2

        assert any(utils._equals_rc(path, contig) for path in paths)
        assert any(utils._equals_rc(path, branch) for path in paths)

    def test_assemble_right_triple_fork(self, right_triple_fork_structure):
        # assemble three contigs from a trip fork
        (graph, contig, L, HDN, R,
         top_sequence, bottom_sequence) = right_triple_fork_structure
        lh = khmer.GraphLabels(graph)
        asm = khmer.SimpleLabeledAssembler(lh)

        hdn = graph.find_high_degree_nodes(contig)
        hdn += graph.find_high_degree_nodes(top_sequence)
        hdn += graph.find_high_degree_nodes(bottom_sequence)
        print(list(hdn))
        lh.label_across_high_degree_nodes(contig, hdn, 1)
        lh.label_across_high_degree_nodes(top_sequence, hdn, 2)
        lh.label_across_high_degree_nodes(bottom_sequence, hdn, 3)
        print(lh.get_tag_labels(list(hdn)[0]))

        paths = asm.assemble(contig[:K])
        print([len(x) for x in paths])

        assert len(paths) == 3

        assert any(utils._equals_rc(path, contig) for path in paths)
        assert any(utils._equals_rc(path, top_sequence) for path in paths)
        assert any(utils._equals_rc(path, bottom_sequence) for path in paths)

    def test_assemble_left_double_fork(self, left_double_fork_structure):
        # assemble entire contig + branch points b/c of labels; start from end
        graph, contig, L, HDN, R, branch = left_double_fork_structure
        lh = khmer.GraphLabels(graph)
        asm = khmer.SimpleLabeledAssembler(lh)

        # first try without the labels
        paths = asm.assemble(contig[-K:])

        assert len(paths) == 1
        # without labels, should get the beginning of the HDN thru the end
        assert paths[0] == contig[HDN.pos:]

        # now add labels and check that we get two full length paths
        hdn = graph.find_high_degree_nodes(contig)
        hdn += graph.find_high_degree_nodes(branch)
        print(list(hdn))
        lh.label_across_high_degree_nodes(contig, hdn, 1)
        lh.label_across_high_degree_nodes(branch, hdn, 2)
        print(lh.get_tag_labels(list(hdn)[0]))

        paths = asm.assemble(contig[-K:])

        assert len(paths) == 2

        assert any(utils._equals_rc(path, contig) for path in paths)
        assert any(utils._equals_rc(path, branch) for path in paths)

    def test_assemble_snp_bubble_single(self, snp_bubble_structure):
        # assemble entire contig + one of two paths through a bubble
        graph, wildtype, mutant, HDN_L, HDN_R = snp_bubble_structure
        lh = khmer.GraphLabels(graph)
        asm = khmer.SimpleLabeledAssembler(lh)

        hdn = graph.find_high_degree_nodes(wildtype)
        assert len(hdn) == 2
        lh.label_across_high_degree_nodes(wildtype, hdn, 1)

        paths = asm.assemble(wildtype[:K])

        assert len(paths) == 1
        assert utils._equals_rc(paths[0], wildtype)

    def test_assemble_snp_bubble_both(self, snp_bubble_structure):
        # assemble entire contig + both paths
        graph, wildtype, mutant, HDN_L, HDN_R = snp_bubble_structure
        lh = khmer.GraphLabels(graph)
        asm = khmer.SimpleLabeledAssembler(lh)

        hdn = graph.find_high_degree_nodes(wildtype)
        hdn += graph.find_high_degree_nodes(mutant)
        assert len(hdn) == 2
        lh.label_across_high_degree_nodes(wildtype, hdn, 1)
        lh.label_across_high_degree_nodes(mutant, hdn, 2)

        paths = asm.assemble(wildtype[:K])

        assert len(paths) == 2

        assert any(utils._contains_rc(wildtype, path) for path in paths)
        assert any(utils._contains_rc(mutant, path) for path in paths)
        # assert all(path[:HDN_L.pos+K][-K:] == HDN_L for path in paths)
        # assert all(path[HDN_R.pos:][:K] == HDN_R for path in paths)
        # assert paths[0][:HDN_L.pos+K] == paths[1][:HDN_L.pos+K]
        # assert paths[0][HDN_R.pos:] == paths[1][HDN_R.pos:]

    def test_assemble_snp_bubble_stopbf(self, snp_bubble_structure):
        # assemble one side of bubble, blocked with stop_filter,
        # when labels on both branches
        # stop_filter should trip a filter failure, negating the label spanning
        graph, wildtype, mutant, HDN_L, HDN_R = snp_bubble_structure
        stop_filter = khmer.Nodegraph(K, 1e5, 4)
        lh = khmer.GraphLabels(graph)
        asm = khmer.SimpleLabeledAssembler(lh, stop_filter=stop_filter)

        hdn = graph.find_high_degree_nodes(wildtype)
        hdn += graph.find_high_degree_nodes(mutant)
        assert len(hdn) == 2
        lh.label_across_high_degree_nodes(wildtype, hdn, 1)
        lh.label_across_high_degree_nodes(mutant, hdn, 2)

        # do the labeling, but block the mutant with stop_filter
        stop_filter.count(mutant[HDN_L.pos + 1:HDN_L.pos + K + 1])
        paths = asm.assemble(wildtype[:K])

        assert len(paths) == 1
        assert any(utils._equals_rc(path, wildtype) for path in paths)

    # @pytest.mark.skip(reason='destroys your computer and then the world')
    def test_assemble_tandem_repeats(self, tandem_repeat_structure):
        # assemble one copy of a tandem repeat
        graph, repeat, tandem_repeats = tandem_repeat_structure
        lh = khmer.GraphLabels(graph)
        asm = khmer.SimpleLabeledAssembler(lh)
        paths = asm.assemble(repeat[:K])

        assert len(paths) == 1
        # There are K-1 k-mers spanning the junction between
        # the beginning and end of the repeat
        assert len(paths[0]) == len(repeat) + K - 1


class TestJunctionCountAssembler:

    def test_beginning_to_end_across_tip(self, right_tip_structure):
        # assemble entire contig, ignoring branch point b/c of labels
        graph, contig, L, HDN, R, tip = right_tip_structure
        asm = khmer.JunctionCountAssembler(graph)
        asm.consume(contig)
        asm.consume(contig)
        asm.consume(contig)

        path = asm.assemble(contig[:K])
        print('P:', path[0])
        print('T:', tip)
        print('C:', contig)
        assert len(path) == 1, "there should only be one path"
        path = path[0]  # @CTB

        assert len(path) == len(contig)
        assert utils._equals_rc(path, contig)
