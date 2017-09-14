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
from khmer import reverse_complement as revcomp
from . import khmer_tst_utils as utils

import pytest
import screed


# We just define this globally rather than in a module-level fixture,
# as we need it during parameterization and whatnot.
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

    To be non-overlapping, the sequence should not include any repeated
    length K-1 k-mers.

    Args:
        exclude (str): If not None, add the k-mers from this sequence to the
        seen set.

    Returns:
        str: A random non-looping sequence.
    '''

    seen = set()

    def add_seen(kmer):
        seen.add(kmer)
        seen.add(revcomp(kmer))

    if exclude is not None:
        for pos in range(0, len(exclude) - K):
            add_seen(exclude[pos:pos + K - 1])

    seq = [random.choice('ACGT') for _ in range(K - 1)]  # do first K-1 bases
    add_seen(''.join(seq))

    while(len(seq) < length):
        next_base = random.choice('ACGT')
        next_kmer = ''.join(seq[-K + 2:] + [next_base])
        assert len(next_kmer) == K - 1
        if (next_kmer) not in seen:
            seq.append(next_base)
            add_seen(next_kmer)
        else:
            continue
    return ''.join(seq)


def reads(sequence, L=100, N=100, dbg_cover=False):
    positions = list(range(len(sequence) - L))
    if dbg_cover is True:
        for start in range(0, len(sequence), K):
            read = sequence[start:start + L]
            if len(read) < K:
                read = sequence[-L:]
            yield read
            N -= 1
    if N < 0:
        return
    for i in range(N):
        start = random.choice(positions)
        yield sequence[start:start + L]


def kmers(sequence):
    for i in range(len(sequence) - K + 1):
        yield sequence[i:i + K]


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


def test_reads():
    contigfile = utils.get_test_data('simple-genome.fa')
    contig = list(screed.open(contigfile))[0].sequence

    for read in reads(contig):
        assert read in contig

    for read in reads(contig):
        assert mutate_sequence(read) not in contig


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
    ~~: Tandem →o→ repeats
'''


@pytest.fixture(params=['simple-genome.fa'])
def known_sequence(request):
    fn = utils.get_test_data(request.param)
    return list(screed.open(fn))[0].sequence


@pytest.fixture(params=list(range(500, 1600, 500)),
                ids=lambda val: '(L={0})'.format(val))
def random_sequence(request):

    def get(exclude=None):
        return get_random_sequence(request.param, exclude=exclude)

    return get


@pytest.fixture(params=[khmer.Nodegraph, khmer.Countgraph],
                ids=['(Type=Nodegraph)', '(Type=Countgraph)'])
def graph(request):

    num_kmers = 50000
    des_fp = 0.00001
    args = optimal_fp(num_kmers, des_fp)
    print('Graph Params:', args)

    return request.param(K, args.htable_size, args.num_htables)


def hdn_counts(sequence, graph):
    '''Get the degree distribution of nodes with degree more than 2.
    '''

    hdns = {}
    for kmer in kmers(sequence):
        d = graph.kmer_degree(kmer)
        if d > 2:
            hdns[d] = hdns.get(d, 0) + 1

    return hdns


@pytest.fixture
def linear_structure(request, graph, random_sequence):
    '''Sets up a simple linear path graph structure.

    sequence
    [0]→o→o~~o→o→[-1]
    '''
    sequence = random_sequence()
    graph.consume(sequence)

    # Check for false positive neighbors in our graph
    # Mark as an expected failure if any are found
    if hdn_counts(sequence, graph):
        request.applymarker(pytest.mark.xfail)

    return graph, sequence


@pytest.fixture(params=[K * 2, -K * 2],
                ids=['(Where={0})'.format(i) for i in ['Start', 'End']])
def right_tip_structure(request, graph, random_sequence):
    '''
    Sets up a graph structure like so:
                                 ([S+1:S+K]+B tip)
    sequence                   ↗
    [0]→o→o~~o→(L)→([S:S+K] HDN)→(R)→o→o→o~~o→[-1]

    Where S is the start position of the high degreen node (HDN).
    That is, it has a single branch at the Sth K-mer.
    '''
    sequence = random_sequence()
    S = request.param
    if S < 0:
        S = len(sequence) + S
    # the HDN
    HDN = Kmer(sequence[S:S + K], pos=S)
    # left of the HDN
    L = Kmer(sequence[S - 1:S - 1 + K], pos=S - 1)
    # right of the HDN
    R = Kmer(sequence[S + 1:S + 1 + K], pos=S + 1)
    # the branch kmer
    tip = Kmer(mutate_position(R, -1),
               pos=R.pos)

    graph.consume(sequence)
    graph.count(tip)

    # Check for false positive neighbors and mark as expected failure if found
    if hdn_counts(sequence, graph) != {3: 1}:
        request.applymarker(pytest.mark.xfail)

    return graph, sequence, L, HDN, R, tip


@pytest.fixture(params=[K * 2, -K * 2],
                ids=['(Where={0})'.format(i) for i in ['Start', 'End']])
def right_double_fork_structure(request, linear_structure, random_sequence):
    '''
    Sets up a graph structure like so:
                                               branch
                                 ([S+1:S+K]+B)→o~~o→o
    core_sequence               ↗
    [0]→o→o~~o→(L)→([S:S+K] HDN)→(R)→o→o→o~~o→[-1]

    Where S is the start position of the high degreen node (HDN)
    and B is the mutated base starting the branch.
    '''

    graph, core_sequence = linear_structure
    print('\nCore Len:', len(core_sequence))
    branch_sequence = random_sequence(exclude=core_sequence)
    print('Branch len:', len(branch_sequence))

    # start position of the HDN
    S = request.param
    if S < 0:
        S = len(core_sequence) + S
    # the HDN
    HDN = Kmer(core_sequence[S:S + K], pos=S)
    # left of the HDN
    L = Kmer(core_sequence[S - 1:S - 1 + K], pos=S - 1)
    # right of the HDN
    R = Kmer(core_sequence[S + 1:S + 1 + K], pos=S + 1)
    # the branch sequence, mutated at position S+1
    branch_start = core_sequence[:R.pos] + mutate_position(R, -1)
    branch_sequence = branch_start + branch_sequence

    graph.consume(core_sequence)
    graph.consume(branch_sequence)

    # Check for false positive neighbors and mark as expected failure if found
    core_hdns = hdn_counts(core_sequence, graph)
    branch_hdns = hdn_counts(branch_sequence, graph)

    # the core and branch sequences should each have exactly
    # ONE node of degree 3 (HDN)
    if core_hdns != {3: 1} or branch_hdns != {3: 1}:
        print(core_hdns, branch_hdns)
        request.applymarker(pytest.mark.xfail)

    return graph, core_sequence, L, HDN, R, branch_sequence


@pytest.fixture
def right_triple_fork_structure(request, right_double_fork_structure,
                                random_sequence):
    '''
    Sets up a graph structure like so:

                                       top_branch
                                ([:S+1]+B)→o~~o→o
    core_sequence              ↗
    [0]→o→o~~o→(L)→([S:S+K] HDN)→(R)→o→o→o~~o→[-1]
                               ↘
                                ([:S+1]+B)→o~~o→o
                                     bottom_branch

    Where S is the start position of the high degreen node (HDN).
    '''

    graph, core_sequence, L, HDN, R, top_sequence = right_double_fork_structure
    bottom_branch = random_sequence(exclude=core_sequence + top_sequence)
    print(len(core_sequence), len(top_sequence), len(bottom_branch))

    # the branch sequence, mutated at position S+1
    # choose a base not already represented at that position
    bases = {'A', 'C', 'G', 'T'}
    mutated = random.choice(list(bases - {R[-1], top_sequence[R.pos + K - 1]}))

    bottom_sequence = core_sequence[:HDN.pos + K] + mutated + bottom_branch

    graph.consume(bottom_sequence)

    # Check for false positive neighbors and mark as expected failure if found
    core_hdns = hdn_counts(core_sequence, graph)
    top_hdns = hdn_counts(top_sequence, graph)
    bottom_hdns = hdn_counts(bottom_sequence, graph)

    # the core, top, and bottom sequences should each have exactly
    # ONE node of degree 4 (HDN)
    if not (core_hdns == top_hdns == bottom_hdns == {4: 1}):
        print(core_hdns, top_hdns, bottom_hdns)
        request.applymarker(pytest.mark.xfail)

    return graph, core_sequence, L, HDN, R, top_sequence, bottom_sequence


@pytest.fixture(params=[K * 2, -K * 2],
                ids=['(Where={0})'.format(i) for i in ['Start', 'End']])
def left_tip_structure(request, graph, random_sequence):
    '''
    Sets up a graph structure like so:

    branch
    (B+[S:S+K-1] tip)
                     ↘                    sequence
        [0]→o~~o→(L)→([S:S+K] HDN)→(R)→o→o~~o→[-1]

    Where S is the start position of the HDN.
    '''
    sequence = random_sequence()
    S = request.param
    if S < 0:
        S = len(sequence) + S
    tip = Kmer(mutate_position(sequence[S - 1:S - 1 + K], 0),
               pos=S - 1 + K)
    HDN = Kmer(sequence[S:S + K], pos=S)
    L = Kmer(sequence[S - 1:S - 1 + K], pos=S - 1)
    R = Kmer(sequence[S + 1:S + 1 + K], pos=S + 1)

    graph.consume(sequence)
    graph.count(tip)

    # Check for false positive neighbors and mark as expected failure if found
    if hdn_counts(sequence, graph) != {3: 1}:
        request.applymarker(pytest.mark.xfail)

    return graph, sequence, L, HDN, R, tip


@pytest.fixture(params=[K * 2, -K * 2],
                ids=['(Where={0})'.format(i) for i in ['Start', 'End']])
def left_double_fork_structure(request, linear_structure, random_sequence):
    '''
    Sets up a graph structure like so:

    o→o~~o→(B+[S:S+K-1])
                        ↘                  core_sequence
          [0]→o→o~~o→(L)→([S:S+K] HDN)→(R)→o→o→o~~o→[-1]

    Where S is the start position of the high degreen node (HDN).
    '''

    graph, core_sequence = linear_structure
    branch_sequence = random_sequence(exclude=core_sequence)

    # start position of the HDN
    S = request.param
    if S < 0:
        S = len(core_sequence) + S
    # the HDN
    HDN = Kmer(core_sequence[S:S + K], pos=S)
    # left of the HDN
    L = Kmer(core_sequence[S - 1:S - 1 + K], pos=S - 1)
    # right of the HDN
    R = Kmer(core_sequence[S + 1:S + 1 + K], pos=S + 1)
    # the branch sequence, mutated at position 0 in L,
    # whih is equivalent to the K-1 prefix of HDN prepended with a new base
    branch_start = mutate_position(L, 0)
    branch_sequence = branch_sequence + \
        branch_start + core_sequence[L.pos + K:]

    graph.consume(core_sequence)
    graph.consume(branch_sequence)

    # Check for false positive neighbors and mark as expected failure if found
    core_hdns = hdn_counts(core_sequence, graph)
    branch_hdns = hdn_counts(branch_sequence, graph)

    # the core and branch sequences should each have exactly
    # ONE node of degree 3 (HDN)
    if not (core_hdns == branch_hdns == {3: 1}):
        request.applymarker(pytest.mark.xfail)

    return graph, core_sequence, L, HDN, R, branch_sequence


@pytest.fixture(params=[K * 2, (-K * 2) - 2],
                ids=['(Where={0})'.format(i) for i in ['Start', 'End']])
def snp_bubble_structure(request, linear_structure):
    '''
    Sets up a graph structure resulting from a SNP (Single Nucleotide
    Polymorphism).

                        (HDN_L[1:]+SNP)→o~~o→(SNP+)
                      ↗                            ↘
    o~~([S:S+K] HDN_L)                             ([S+K+1:S+2K+1] HDN_R)~~o
                      ↘                           ↗
                        (HDN_L[1:]+W)→o~~o~~o→(W+)

    Where S is the start position of HDN directly left of the SNP (HDN_L),
    SNP is the mutated base, and W is the wildtype (original) base.
    Of course, W and SNP could be interchanged here, we don't actually
    know which is which ;)

    Note our parameterization: we need a bit more room from the ends,
    so we bring the rightmost SNP a tad left.
    '''

    graph, wildtype_sequence = linear_structure
    S = request.param
    if S < 0:
        S = len(wildtype_sequence) + S
    snp_sequence = mutate_position(wildtype_sequence, S + K)
    HDN_L = Kmer(wildtype_sequence[S:S + K], pos=S)
    HDN_R = Kmer(wildtype_sequence[S + K + 1:S + 2 * K + 1], pos=S + K + 1)

    graph.consume(wildtype_sequence)
    graph.consume(snp_sequence)

    # Check for false positive neighbors and mark as expected failure if found
    w_hdns = hdn_counts(wildtype_sequence, graph)
    snp_hdns = hdn_counts(snp_sequence, graph)
    if not (w_hdns == snp_hdns == {3: 2}):
        print(w_hdns, snp_hdns)
        print(HDN_L, HDN_R)
        print(wildtype_sequence[HDN_L.pos + K + 1])
        print(snp_sequence[HDN_L.pos + K + 1])
        request.applymarker(pytest.mark.xfail)

    return graph, wildtype_sequence, snp_sequence, HDN_L, HDN_R


@pytest.fixture(params=[2, 3, 4, 5, 6, 7, 8])
def tandem_repeat_structure(request, linear_structure):

    graph, sequence = linear_structure

    tandem_repeats = sequence * request.param
    graph.consume(tandem_repeats)

    if hdn_counts(tandem_repeats, graph):
        request.applymarker(pytest.mark.xfail)

    return graph, sequence, tandem_repeats


@pytest.fixture
def circular_linear_structure(request, linear_structure):
    graph, sequence = linear_structure

    sequence += sequence

    if hdn_counts(sequence, graph):
        request.applymarker(pytest.mark.xfail)

    return graph, sequence
