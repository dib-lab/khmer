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

K = 21

class Kmer(str):

    def __init__(self, value, pos=0):
        self.pos = pos

    def __new__(cls, value, pos=0):
        if not len(value) == K:
            raise ValueError('bad k-mer length')
        return str.__new__(cls, value)


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


def get_random_sequence(length, exclude=None):
    '''Generate a random (non-looping) nucleotide sequence.

    Args:
        exlcude (str): If not None, add the k-mers from this sequence to the seen set.
    
    Returns:
        str: A random non-looping sequence.
    '''
    seen = set()
    if exclude is not None:
        print('Adding', exclude, 'to seen set')
        for pos in range(0, len(exclude)-K+1):
            seen.add(exclude[pos:pos+K])

    seq = [random.choice('ACGT') for _ in range(K)] # do first K bases
    seen.add(''.join(seq))
    while(len(seq) < length):
        next_base = random.choice('ACGT')
        next_kmer = ''.join(seq[-K+1:] + [next_base])
        assert len(next_kmer) == K
        if (next_kmer) not in seen:
            seq.append(next_base)
        else:
            continue
    return ''.join(seq)  


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


@pytest.fixture(params=['simple-genome.fa'])
def known_sequence(request):
    fn = utils.get_test_data(request.param)
    return list(screed.open(fn))[0].sequence


@pytest.fixture(params=list(range(500, 1500, 100)))
def random_sequence(request):

    return get_random_sequence(request.param) 


@pytest.fixture(params=[khmer.Nodegraph, khmer.Countgraph])
def graph(request):
    return request.param(K, 1e6, 4)



'''
# GRAPH STRUCTURE FIXTURES

These fixtures emit various graph structures with their corresponding
sequences and important nodes. They take a random sequence fixture and
a graph fixture, then consume sequence and generate k-mers accordingly.

We're using a bespoke but simple language to describe graph structures in the
docstrings of these tests. It is as follows:

    o: Node
    [x:y]: Node at position in sequence
    [x:y]+S: Node at position in sequence with extra base (where S in ACGT)
    (Name), ([x:y] Name): Named node, named node at position
    → : Edge
    ~~: Tandem → o→  repeats

'''

@pytest.fixture
def linear_structure(graph, random_sequence):
    '''Sets up a simple linear path graph structure.

    contig
    [0]→ o→ o~~o→ o→ [-1]
    '''

    graph.consume(random_sequence)

    return graph, random_sequence


@pytest.fixture(params=[K*2, K*3, -K*3, -K*2])
def right_tip_structure(request, graph, random_sequence):
    '''
    Sets up a graph structure like so:
                                        ([S+1:S+K]+B tip)
    contig                             ↗
    [0]→ o→ o~~o→ (L)→ ([S:S+K] HDN)→ (R)→ o→ o→ o~~o→ [-1]

    Where S is the start position of the high degreen node (HDN).
    That is, it has a single branch at the Sth K-mer.
    '''
    
    S = request.param
    if S < 0:
        S = len(random_sequence) + S
    # the HDN
    HDN = Kmer(random_sequence[S:S+K], pos=S)
    # the branch kmer
    tip = Kmer(mutate_position(random_sequence[S+1:S+1+K], -1),
               pos=S+1)
    # left of the HDN
    L = Kmer(random_sequence[S-1:S-1+K], pos=S-1)
    # right of the HDN
    R = Kmer(random_sequence[S+1:S+1+K], pos=S+1)


    graph.consume(random_sequence)
    graph.count(tip)

    return graph, random_sequence, L, HDN, R, tip


@pytest.fixture(params=[K*2, K*3, -K*3, -K*2])
def right_fork_structure(request, graph, random_sequence):
    '''
    Sets up a graph structure like so:
                                        ([:S+1]+B*25 branch)
    contig                             ↗
    [0]→ o→ o~~o→ (L)→ ([S:S+K] HDN)→ (R)→ o→ o→ o~~o→ [-1]

    Where S is the start position of the high degreen node (HDN).
    The branch is fixed at length 25, and will not form a loop
    with the main contig.
    '''
    
    S = request.param
    if S < 0:
        S = len(random_sequence) + S
    # the HDN
    HDN = Kmer(random_sequence[S:S+K], pos=S)
    # the branch sequence, mutated at position S+1
    branch = mutate_position(random_sequence[:S+2], -1)
    branch += get_random_sequence(25, exclude=random_sequence)
    # left of the HDN
    L = Kmer(random_sequence[S-1:S-1+K], pos=S-1)
    # right of the HDN
    R = Kmer(random_sequence[S+1:S+1+K], pos=S+1)


    graph.consume(random_sequence)
    graph.consume(branch)

    return graph, random_sequence, L, HDN, R, branch


@pytest.fixture(params=[K*2, K*3, -K*3, -K*2])
def left_tip_structure(request, graph, random_sequence):
    '''
    Sets up a graph structure like so:

    (B+[S:S+K-1] tip)
                     ↘
    [0]→ o~~o→ (L)→ ([S:S+K] HDN)→ (R)→ o→ o~~o→ [-1]

    Where S is the start position of the HDN
    That is, it has a single branch at the 100th K-mer.
    '''
    S = request.param
    if S < 0:
        S = len(random_sequence) + S
    tip = Kmer(mutate_position(random_sequence[S-1:S-1+K], 0),
               pos=S-1+K)
    HDN = Kmer(random_sequence[S:S+K], pos=S)
    L = Kmer(random_sequence[S-1:S-1+K], pos=S-1)
    R = Kmer(random_sequence[S+1:S+1+K], pos=S+1)

    graph.consume(random_sequence)
    graph.count(tip)

    return graph, random_sequence, L, HDN, R, tip


class TestNonBranching:

    def test_all_start_positions(self, linear_structure):
        # assemble entire contig, starting from wherever
        nodegraph, contig = linear_structure

        for start in range(0, len(contig), 150):
            path = nodegraph.assemble_linear_path(contig[start:start + K])
            assert utils._equals_rc(path, contig), start


class TestLinearAssembler_RightBranching:
 
    def test_branch_point(self, right_tip_structure):
        graph, contig, L, HDN, R, tip = right_tip_structure

        assert graph.kmer_degree(HDN) == 3

    def test_beginning_to_branch(self, right_tip_structure):
        # assemble from beginning of contig, up until branch point
        graph, contig, L, HDN, R, tip = right_tip_structure

        path = graph.assemble_linear_path(contig[0:K])

        assert len(path) == HDN.pos + K
        assert utils._equals_rc(path, contig[:len(path)])

    def test_beginning_to_branch_revcomp(self, right_tip_structure):
        # assemble from beginning of contig, up until branch point
        # starting from rev comp
        graph, contig, L, HDN, R, tip = right_tip_structure
        path = graph.assemble_linear_path(revcomp(contig[0:K]))

        assert len(path) == HDN.pos + K
        assert utils._equals_rc(path, contig[:len(path)])

    def test_left_of_branch_to_beginning(self, right_tip_structure):
        # start from HDN (left of branch)
        graph, contig, L, HDN, R, tip = right_tip_structure
        path = graph.assemble_linear_path(L)

        assert len(path) == HDN.pos+K
        assert utils._equals_rc(path, contig[:len(path)])

    def test_left_of_branch_to_beginning_revcomp(self, right_tip_structure):
        # start from revcomp of HDN (left of branch)
        graph, contig, L, HDN, R, tip = right_tip_structure
        path = graph.assemble_linear_path(revcomp(L))

        assert len(path) == HDN.pos+K
        assert utils._equals_rc(path, contig[:len(path)])

    def test_right_of_branch_outwards_to_ends(self, right_tip_structure):
        # assemble from right of branch point (at R)
        # Should get the *entire* original contig, as the assembler
        # will move left relative to the branch, and not consider it
        # as a high degree node
        graph, contig, L, HDN, R, tip = right_tip_structure
        path = graph.assemble_linear_path(R)

        assert len(path) == len(contig)
        assert utils._equals_rc(path, contig)

    def test_end_to_beginning(self, right_tip_structure):
        # should have exact same behavior as right_of_branch_outwards
        graph, contig, L, HDN, R, tip = right_tip_structure
        path = graph.assemble_linear_path(contig[-K:])

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

        path = graph.assemble_linear_path(contig[-K:])

        assert len(path) == len(contig) - HDN.pos
        assert utils._equals_rc(path, contig[HDN.pos:])

    def test_branch_to_end(self, left_tip_structure):
        # assemble from branch point until end
        graph, contig, L, HDN, R, tip = left_tip_structure

        path = graph.assemble_linear_path(HDN)

        assert len(path) == len(contig) - HDN.pos
        assert utils._equals_rc(path, contig[HDN.pos:])

    def test_branch_outwards_to_ends_with_stopbf(self, left_tip_structure):
        # block the tip with the stop_bf. should return a full length contig.
        graph, contig, L, HDN, R, tip = left_tip_structure

        stop_bf = khmer.Nodegraph(K, 1e5, 4)
        stop_bf.count(tip) 

        path = graph.assemble_linear_path(HDN, stop_bf)

        assert len(path) == len(contig)
        assert utils._equals_rc(path, contig)

    def test_branch_outwards_to_ends_with_stopbf_revcomp(self, left_tip_structure):
         # block the tip with the stop_bf. should return a full length contig.
        graph, contig, L, HDN, R, tip = left_tip_structure

        stop_bf = khmer.Nodegraph(K, 1e5, 4)
        stop_bf.count(tip) 

        path = graph.assemble_linear_path(revcomp(HDN), stop_bf)

        assert len(path) == len(contig)
        assert utils._equals_rc(path, contig)


    def test_end_thru_tip_with_stopbf(self, left_tip_structure):
        # assemble up to branch point, and include introduced branch b/c
        # of stop bf
        graph, contig, L, HDN, R, tip = left_tip_structure

        stop_bf = khmer.Nodegraph(K, 1e5, 4)
        stop_bf.count(L)          # ...and block original path
        path = graph.assemble_linear_path(contig[-K:], stop_bf)

        assert len(path) == len(contig) - HDN.pos + 1

        # should be the tip k-kmer, plus the last base of the HDN thru
        # the end of the contig
        assert utils._equals_rc(path, tip + contig[HDN.pos+K-1:])


    def test_single_node_flanked_by_hdns(self, left_tip_structure):
        # assemble single node flanked by high-degree nodes
        # we'll copy the main nodegraph before mutating it
        graph, contig, L, HDN, R, tip = left_tip_structure

        graph.consume(mutate_position(contig, HDN.pos + K))

        path = graph.assemble_linear_path(HDN)

        assert len(path) == K
        assert utils._equals_rc(path, HDN)


class TestLabeledAssembler:

    def test_beginning_to_end_across_tip(self, right_tip_structure):
        # assemble entire contig, ignoring branch point b/c of labels
        graph, contig, L, HDN, R, tip = right_tip_structure
        lh = khmer._GraphLabels(graph)

        hdn = graph.find_high_degree_nodes(contig)
        # L, HDN, and R will be labeled with 1
        lh.label_across_high_degree_nodes(contig, hdn, 1)

        path = lh.assemble_labeled_path(contig[:K])
        assert len(path) == 1, "there should only be one path"
        path = path[0]  # @CTB

        assert len(path) == len(contig)
        assert utils._equals_rc(path, contig)


    def test_assemble_original_and_branch_contigs(self, right_fork_structure):
        # assemble entire contig + branch point b/c of labels
        graph, contig, L, HDN, R, branch = right_fork_structure
        lh = khmer._GraphLabels(graph)

        hdn = graph.find_high_degree_nodes(contig)
        hdn += graph.find_high_degree_nodes(branch)
        print(list(hdn))
        lh.label_across_high_degree_nodes(contig, hdn, 1)
        lh.label_across_high_degree_nodes(branch, hdn, 2)
        print(lh.get_tag_labels(list(hdn)[0]))

        paths = lh.assemble_labeled_path(contig[:K])
        print('Path lengths', [len(x) for x in paths])

        assert len(paths) == 2

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
