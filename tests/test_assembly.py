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
from khmer import ReadParser
from khmer import reverse_complement as revcomp

import screed

import pytest

from . import khmer_tst_utils as utils


def teardown():
    utils.cleanup()

def setup_module(module):
    module.K = 21

def mutate_base(base):
    if base in 'AT':
        return random.choice('GC')
    elif base in 'GC':
        return random.choice('AT')
    else:
        assert False, 'bad base'
    
def mutate_sequence(sequence, N=1):
    sequence = list(sequence)
    positions = random.sample(range(len(sequence)), N)
    
    for i in positions:
        sequence[i] = mutate_base(sequence[i])
        
    return ''.join(sequence)

def mutate_position(sequence, pos):
    sequence = list(sequence)
    sequence[pos] = mutate_base(sequence[pos])
    return ''.join(sequence)


def reads_from_sequence(sequence, L=100, N=100):
    positions = list(range(len(sequence) - L))
    for i in range(N):
        start = random.choice(positions)
        yield sequence[start:start+L]


def test_mutate_sequence():
    for _ in range(100):
        assert 'A' not in mutate_sequence('A' * 10, 10)
        assert 'T' not in mutate_sequence('T' * 10, 10)
        assert 'C' not in mutate_sequence('C' * 10, 10)
        assert 'G' not in mutate_sequence('G' * 10, 10)


def test_mutate_position():
    assert mutate_position('AAAA', 2) in ['AACA', 'AAGA']
    assert mutate_position('TTTT', 2) in ['TTCT', 'TTGT']
    assert mutate_position('CCCC', 2) in ['CCAC', 'CCTC']
    assert mutate_position('GGGG', 2) in ['GGAG', 'GGTG']


def test_reads_from_sequence():
    contigfile = utils.get_test_data('simple-genome.fa')
    contig = list(screed.open(contigfile))[0].sequence

    for read in reads_from_sequence(contig):
        assert read in contig
        
    for read in reads_from_sequence(contig):
        assert mutate_sequence(read) not in contig


class TestNonBranching:
    '''
    Sets up a simple linear path graph structure.

    contig
    [0]--o--o--o~~o--o--o--[-1]
    '''

    def setup_class(self):
        contigfile = utils.get_test_data('simple-genome.fa')
        self.contig = list(screed.open(contigfile))[0].sequence

    def test_all_start_positions(self):
        # assemble entire contig, starting from wherever
        print('contig len', len(self.contig))

        nodegraph = khmer.Nodegraph(K, 1e5, 4)

        nodegraph.consume(self.contig)

        for start in range(0, len(self.contig), 150):
            path = nodegraph.assemble_linear_path(self.contig[start:start + K])
            assert utils._equals_rc(path, self.contig), start


class TestLinearAssembler_RightBranching:
    '''
    Sets up a graph structure like so:
                                 [101:101+K-1]+G(tip)
    contig                     ↗
    [0]→ o→ o→ o→ L→ [100:100+K]→ R→ o→ o→ o~~o→ [-1]

    That is, it has a single branch at the 100th K-mer.
    '''

    def setup_class(self):
        contigfile = utils.get_test_data('simple-genome.fa')

        self.contig = list(screed.open(contigfile))[0].sequence
        self.tip = self.contig[101:101+K-1] + 'G'
        self.HDN = self.contig[100:100+K]

        self.nodegraph = khmer.Nodegraph(K, 1e5, 4)
        self.nodegraph.consume(self.contig)
        self.nodegraph.count(self.tip)

    def test_branch_point(self):
        assert self.nodegraph.kmer_degree(self.HDN) == 3

    def test_beginning_to_branch(self):
        # assemble from beginning of contig, up until branch point

        path = self.nodegraph.assemble_linear_path(self.contig[0:K])
        len_path = len(path)

        assert len_path == 100+K
        assert utils._equals_rc(path, self.contig[:len_path])

    def test_beginning_to_branch_revcomp(self):
        # assemble from beginning of contig, up until branch point
        # starting from rev comp
        path = self.nodegraph.assemble_linear_path(revcomp(self.contig[0:K]))
        len_path = len(path)

        assert utils._equals_rc(path, self.contig[:len_path])

    def test_left_of_branch_to_beginning(self):
        # start from HDN (left of branch)
        path = self.nodegraph.assemble_linear_path(self.HDN)
        len_path = len(path)

        assert len_path == 100+K
        assert utils._equals_rc(path, self.contig[:len_path])

    def test_left_of_branch_to_beginning_revcomp(self):
        # start from revcomp of HDN (left of branch)
        path = self.nodegraph.assemble_linear_path(self.HDN)
        len_path = len(path)

        assert len_path == 100+K
        assert utils._equals_rc(path, self.contig[:len_path])

    def test_right_of_branch_outwards_to_ends(self):
        # assemble from right of branch point (at R)
        # Should get the *entire* original contig, as the assembler
        # will move left relative to the branch, and not consider it
        # as a high degree node
        path = self.nodegraph.assemble_linear_path(self.contig[101:101 + K])
        len_path = len(path)

        assert len_path == len(self.contig)
        assert utils._equals_rc(path, self.contig)

    def test_end_to_beginning(self):
        # should have exact same behavior as right_of_branch_outwards
        path = self.nodegraph.assemble_linear_path(self.contig[-K:])
        len_path = len(path)

        assert len_path == len(self.contig)
        assert utils._equals_rc(path, self.contig)


class TestLinearAssembler_LeftBranching:
    '''
    Sets up a graph structure like so:

    T+[101:101+K-1](tip)
                       ↘
    [0]→ o→ o→ o→ L→ [101:101+K]→ R→ o→ o→ o~~o→ [-1]

    That is, it has a single branch at the 100th K-mer.
    '''

    def setup_class(self):
        contigfile = utils.get_test_data('simple-genome.fa')

        self.contig = list(screed.open(contigfile))[0].sequence
        self.tip = 'T' + self.contig[101:101+K-1]
        self.HDN = self.contig[101:101+K]
        self.L = self.contig[100:100+K]
        self.R = self.contig[102:102+K]

        self.nodegraph = khmer.Nodegraph(K, 1e5, 4)
        self.nodegraph.consume(self.contig)
        self.nodegraph.count(self.tip)

    def test_branch_point(self):
        assert self.nodegraph.kmer_degree(self.HDN) == 3

    def test_end_to_branch(self):
        # assemble from end until branch point
        # should include HDN
        path = self.nodegraph.assemble_linear_path(self.contig[-K:])
        len_path = len(path)

        assert len_path == 899
        assert utils._equals_rc(path, self.contig[101:])

    def test_branch_to_end(self):
        # assemble from branch point until end
        path = self.nodegraph.assemble_linear_path(self.HDN)
        len_path = len(path)

        assert len_path == 899
        assert utils._equals_rc(path, self.contig[101:])

    def test_branch_outwards_to_ends_with_stopbf(self):
        # block the tip with the stop_bf. should return a full length contig.
        stop_bf = khmer.Nodegraph(K, 1e5, 4)
        stop_bf.count(self.tip) 

        path = self.nodegraph.assemble_linear_path(self.HDN, stop_bf)
        len_path = len(path)

        assert len_path == 1000
        assert utils._equals_rc(path, self.contig)

    def test_branch_outwards_to_ends_with_stopbf_revcomp(self):
         # block the tip with the stop_bf. should return a full length contig.
        stop_bf = khmer.Nodegraph(K, 1e5, 4)
        stop_bf.count(self.tip) 

        path = self.nodegraph.assemble_linear_path(revcomp(self.HDN), stop_bf)
        len_path = len(path)

        assert len_path == 1000
        assert utils._equals_rc(path, self.contig)


    def test_end_thru_tip_with_stopbf(self):
        # assemble up to branch point, and include introduced branch b/c
        # of stop bf
 
        stop_bf = khmer.Nodegraph(K, 1e5, 4)
        stop_bf.count(self.L)          # ...and block original path

        path = self.nodegraph.assemble_linear_path(self.contig[-K:], stop_bf)
        len_path = len(path)

        assert len_path == 900
        assert utils._equals_rc(path, 'T' + self.contig[101:])


    def test_single_node_flanked_by_hdns(self):
        # assemble single node flanked by high-degree nodes
        # we'll copy the main nodegraph before mutating it
        nodegraph = khmer.Nodegraph(K, 1e5, 4)
        nodegraph.update(self.nodegraph)

        # add a neighbor after the current tip
        nodegraph.count(self.contig[102:102+K-1] + 'G')

        path = nodegraph.assemble_linear_path(self.HDN)

        assert len(path) == K
        assert utils._equals_rc(path, self.HDN)


def test_assemble_linear_path_single_node_interrupted():
    # assemble single node.
    contigfile = utils.get_test_data('simple-genome.fa')
    contig = list(screed.open(contigfile))[0].sequence
    print('contig len', len(contig))

    K = 21

    nodegraph = khmer.Nodegraph(K, 1e5, 4)

    nodegraph.consume(contig)
    nodegraph.consume(contig[:110] + 'G')  # will add another neighbor/middle

    path = nodegraph.assemble_linear_path(contig[100:121])
    len_path = len(path)

    print('len path:', len_path)

    assert utils._equals_rc(path, contig)       # this is bad behavior...


def test_assemble_labeled_paths():
    # assemble entire contig, ignoring branch point b/c of labels
    contigfile = utils.get_test_data('simple-genome.fa')
    contig = list(screed.open(contigfile))[0].sequence
    print('contig len', len(contig))

    K = 21

    nodegraph = khmer.Nodegraph(K, 1e5, 4)
    lh = khmer._GraphLabels(nodegraph)

    nodegraph.consume(contig)
    nodegraph.count(contig[100:120] + 'T')  # will add another neighbor

    print(contig[100:125])

    hdn = nodegraph.find_high_degree_nodes(contig)
    lh.label_across_high_degree_nodes(contig, hdn, 1)

    path = lh.assemble_labeled_path(contig[:K])
    path = path[0]  # @CTB
    len_path = len(path)

    print('len path:', len_path)

    assert utils._equals_rc(path, contig)


def test_assemble_labeled_paths_2():
    # assemble entire contig + branch point b/c of labels
    contigfile = utils.get_test_data('simple-genome.fa')
    contig = list(screed.open(contigfile))[0].sequence
    print('contig len', len(contig))

    K = 21

    nodegraph = khmer.Nodegraph(K, 1e5, 4)
    lh = khmer._GraphLabels(nodegraph)

    nodegraph.consume(contig)
    branch = contig[:120] + 'TGATGGACAG'
    nodegraph.consume(branch)  # will add a branch

    hdn = nodegraph.find_high_degree_nodes(contig)
    hdn += nodegraph.find_high_degree_nodes(branch)
    print(list(hdn))
    lh.label_across_high_degree_nodes(contig, hdn, 1)
    lh.label_across_high_degree_nodes(branch, hdn, 2)
    print(lh.get_tag_labels(list(hdn)[0]))

    paths = lh.assemble_labeled_path(contig[:K])
    print([len(x) for x in paths])
    len_path = len(paths)

    print('len path:', len_path)

    found = False
    for path in paths:
        if utils._equals_rc(path, contig):
            found = True
            break
    assert found

    found = False
    for path in paths:
        if utils._equals_rc(path, branch):
            found = True
            break
    assert found


def test_assemble_labeled_paths_3():
    # assemble entire contig + branch points b/c of labels
    contigfile = utils.get_test_data('simple-genome.fa')
    contig = list(screed.open(contigfile))[0].sequence
    print('contig len', len(contig))

    K = 21

    nodegraph = khmer.Nodegraph(K, 1e5, 4)
    lh = khmer._GraphLabels(nodegraph)

    nodegraph.consume(contig)
    branch = contig[:120] + 'TGATGGACAG'
    nodegraph.consume(branch)  # will add a branch
    branch2 = contig[:120] + 'GCGGATGGATGGAGCCGAT'
    nodegraph.consume(branch2)  # will add a third branch

    hdn = nodegraph.find_high_degree_nodes(contig)
    hdn += nodegraph.find_high_degree_nodes(branch)
    hdn += nodegraph.find_high_degree_nodes(branch2)
    print(list(hdn))
    lh.label_across_high_degree_nodes(contig, hdn, 1)
    lh.label_across_high_degree_nodes(branch, hdn, 2)
    lh.label_across_high_degree_nodes(branch2, hdn, 3)
    print(lh.get_tag_labels(list(hdn)[0]))

    paths = lh.assemble_labeled_path(contig[:K])
    print([len(x) for x in paths])
    len_path = len(paths)

    print('len path:', len_path)

    found = False
    for path in paths:
        if utils._equals_rc(path, contig):
            found = True
            break
    assert found

    found = False
    for path in paths:
        if utils._equals_rc(path, branch):
            found = True
            break
    assert found

    found = False
    for path in paths:
        if utils._equals_rc(path, branch2):
            found = True
            break
    assert found


def test_assemble_labeled_paths_4():
    # assemble entire contig + branch points b/c of labels; start from end
    contigfile = utils.get_test_data('simple-genome.fa')
    contig = list(screed.open(contigfile))[0].sequence
    print('contig len', len(contig))

    K = 21

    nodegraph = khmer.Nodegraph(K, 1e5, 4)
    lh = khmer._GraphLabels(nodegraph)

    nodegraph.consume(contig)
    branch = 'TGATGGACAG' + contig[120:]
    nodegraph.consume(branch)  # will add a branch
    branch2 = 'GCGGATGGATGGAGCCGAT' + contig[120:]
    nodegraph.consume(branch2)  # will add a third branch

    hdn = nodegraph.find_high_degree_nodes(contig)
    hdn += nodegraph.find_high_degree_nodes(branch)
    hdn += nodegraph.find_high_degree_nodes(branch2)
    print(list(hdn))
    lh.label_across_high_degree_nodes(contig, hdn, 1)
    lh.label_across_high_degree_nodes(branch, hdn, 2)
    lh.label_across_high_degree_nodes(branch2, hdn, 3)
    print(lh.get_tag_labels(list(hdn)[0]))

    paths = lh.assemble_labeled_path(contig[-K:])
    print([len(x) for x in paths])
    len_path = len(paths)

    print('len path:', len_path)

    found = False
    for path in paths:
        if utils._equals_rc(path, contig):
            found = True
            break
    assert found

    found = False
    for path in paths:
        if utils._equals_rc(path, branch):
            found = True
            break
    assert found

    found = False
    for path in paths:
        if utils._equals_rc(path, branch2):
            found = True
            break
    assert found


def test_assemble_labeled_paths_5():
    # assemble entire contig + one of two paths through a bubble
    contigfile = utils.get_test_data('simple-genome.fa')
    contig = list(screed.open(contigfile))[0].sequence
    print('contig len', len(contig))

    K = 21

    nodegraph = khmer.Nodegraph(K, 1e5, 4)
    lh = khmer._GraphLabels(nodegraph)

    nodegraph.consume(contig)
    contig2 = contig[:200] + 'G' + contig[201:]
    nodegraph.consume(contig2)

    hdn = nodegraph.find_high_degree_nodes(contig)
    assert len(hdn) == 2
    lh.label_across_high_degree_nodes(contig, hdn, 1)

    path = lh.assemble_labeled_path(contig[:K])
    path = path[0]  # @CTB
    len_path = len(path)

    print('len path:', len_path)

    assert utils._equals_rc(path, contig)
