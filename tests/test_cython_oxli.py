from __future__ import print_function
from __future__ import absolute_import

import gc
import itertools
import random

import khmer
from khmer import _oxli
from khmer.khmer_args import estimate_optimal_with_K_and_f as optimal_fp
from khmer import reverse_complement as revcomp
from khmer import reverse_hash as revhash
from . import khmer_tst_utils as utils
from .graph_features import *

import pytest
import screed


def teardown():
    utils.cleanup()


@pytest.fixture
def partitioner(graph):
    sp = _oxli.StreamingPartitioner(graph)
    return graph, sp


@pytest.fixture
def single_component(partitioner, random_sequence):
    graph, partitioner = partitioner
    sequence = random_sequence()
    partitioner.consume_sequence(sequence)
    return graph, partitioner, sequence


class TestStreamingPartitionerBasic:

    def teardown_method(self, method):
        # Force garbage to collect. When Python component objects exist and
        # their underlying c++ Component objects are destroyed, the Python
        # wrapper becomes the sole owner of the pointer. By manually collecting
        # garbage between tests we assure that these objects are freed, and we
        # can properly test the _n_destroyed property to make sure there are no
        # real memory leaks.
        gc.collect()

    def test_one_component(self, known_sequence):
        inpath = utils.get_test_data('random-20-a.fa')

        cg = khmer.Countgraph(K, 1e5, 4)
        sp = _oxli.StreamingPartitioner(cg)
        sp.consume_sequence(known_sequence)

        assert sp.n_components == 1

    def test_two_components(self, random_sequence):
        comp1 = random_sequence()
        comp2 = random_sequence(exclude=comp1)

        cg = khmer.Nodegraph(K, 1e5, 4)
        sp = _oxli.StreamingPartitioner(cg)
        
        sp.consume_sequence(comp1)
        assert sp.n_components == 1
        
        sp.consume_sequence(comp2)
        assert sp.n_components == 2

    def test_components_iter(self, random_sequence):
        comp1 = random_sequence()
        comp2 = random_sequence(exclude=comp1)

        cg = khmer.Nodegraph(K, 1e5, 4)
        sp = _oxli.StreamingPartitioner(cg)
        
        sp.consume_sequence(comp1)
        sp.consume_sequence(comp2)
        assert sp.n_components == 2

        comps = list(sp.components())
        assert len(comps) == 2

    def test_component_n_tags(self, random_sequence):
        seq = random_sequence()

        cg = khmer.Nodegraph(K, 1e5, 4)
        sp = _oxli.StreamingPartitioner(cg)
        sp.consume_sequence(seq)

        tags = [t for t,c in sp.tag_components()]
        comp = sp.get_nearest_component(seq[:K])
        assert len(tags) == len(comp)

    def test_tag_components_iter(self, random_sequence):
        comp1 = random_sequence()
        comp2 = random_sequence(exclude=comp1)

        cg = khmer.Nodegraph(K, 1e5, 4)
        sp = _oxli.StreamingPartitioner(cg)
        
        sp.consume_sequence(comp1)
        sp.consume_sequence(comp2)
        assert sp.n_components == 2

        tags = []
        comps = set()
        for tag, comp in sp.tag_components():
            tags.append(tag)
            comps.add(comp)
        
        assert sum([len([tag for tag in comp]) for comp in comps]) == len(tags)
        assert len(comps) == 2
        assert len(tags) == sum([len(c) for c in comps])

    def test_get_nearest_component(self, random_sequence):
        seq1 =  random_sequence()
        seq2 = random_sequence(exclude=seq1)

        cg = khmer.Nodegraph(K, 1e5, 4)
        sp = _oxli.StreamingPartitioner(cg)
        
        sp.consume_sequence(seq1)
        sp.consume_sequence(seq2)

        c1 = sp.get_nearest_component(seq1[:K])
        c2 = sp.get_nearest_component(seq2[:K])
        assert c1.component_id != c2.component_id

        for tag in c1:
            assert utils._contains_rc(seq1, revhash(tag, K))
            assert not utils._contains_rc(seq2, revhash(tag, K))

        for tag in c2:
            assert utils._contains_rc(seq2, revhash(tag, K))
            assert not utils._contains_rc(seq1, revhash(tag, K))

    def test_merge_components(self, random_sequence):
        seq1 = random_sequence()
        seq2 = random_sequence(exclude=seq1)

        cg = khmer.Nodegraph(K, 1e5, 4)
        sp = _oxli.StreamingPartitioner(cg)
        
        sp.consume_sequence(seq1)
        sp.consume_sequence(seq2)
        assert sp.n_components == 2

        sp.consume_sequence(seq1 + seq2)
        assert sp.n_components == 1

        comps = list(sp.components())
        assert len(comps) == 1


    def test_multi_merge_components(self, random_sequence):
        seq1 = random_sequence()
        seq2 = random_sequence(exclude=seq1)
        seq3 = random_sequence(exclude=seq1+seq2)

        cg = khmer.Nodegraph(K, 1e5, 4)
        sp = _oxli.StreamingPartitioner(cg)
        
        sp.consume_sequence(seq1)
        sp.consume_sequence(seq2)
        sp.consume_sequence(seq3)
        assert sp.n_components == 3

        sp.consume_sequence(seq1 + seq2 + seq3)
        assert sp.n_components == 1

    def test_nomerge_k_minus_2_overlap(self, single_component, random_sequence):
        '''Test that components are not merged when they have a length K-2 overlap.
        '''

        graph, partitioner, seq = single_component
        first = seq[:K-2]
        neighbor = random_sequence(exclude=seq) + first

        assert partitioner.n_components == 1
        partitioner.consume_sequence(neighbor)
        print(seq, neighbor, graph.assemble_linear_path(seq[:K]), sep='\n')
        assert partitioner.n_components == 2

    @pytest.mark.parametrize("where", ["beginning", "end"])
    def test_merge_k_minus_1_overlap(self, single_component, random_sequence,
                                     where):
        '''Test that components are merged when they have a length K-1 overlap.
        '''

        graph, partitioner, seq = single_component
        if where == "beginning":
            overlap = seq[:K-1]
            neighbor = random_sequence(exclude=seq) + overlap
        else:
            overlap = seq[-K+1:]
            neighbor = overlap + random_sequence(exclude=seq)

        assert partitioner.n_components == 1
        partitioner.consume_sequence(neighbor)
        path = graph.assemble_linear_path(seq[:K])
        assert partitioner.n_components == 1

    def test_merge_k_overlap(self, single_component, random_sequence):
        '''Test that components are merged when they have a length K overlap.
        '''

        graph, partitioner, seq = single_component
        first = seq[:K]
        neighbor = random_sequence(exclude=seq) + first

        assert partitioner.n_components == 1
        partitioner.consume_sequence(neighbor)
        print(seq, neighbor, graph.assemble_linear_path(seq[:K]), sep='\n')
        assert partitioner.n_components == 1
        

    @pytest.mark.parametrize("n_reads", list(range(100, 1001, 100)))
    def test_one_component_from_reads(self, random_sequence, n_reads):
        seq = random_sequence()
        seq_reads = list(reads(seq, dbg_cover=True, N=n_reads))

        G = khmer.Nodegraph(K, 1e6, 4)
        sp = _oxli.StreamingPartitioner(G)
        for read in seq_reads:
            sp.consume_sequence(read)

        assert sp.n_components == 1

    @pytest.mark.parametrize("n_components", list(range(1, 10)))
    def test_streaming_multicomponents(self, random_sequence, n_components):
        '''Test with many components from reads, and check for memory leaks.'''
        seqs = []
        for _ in range(n_components):
            seqs.append(random_sequence(exclude=''.join(seqs)))

        seq_reads = []
        for seq in seqs:
            seq_reads.extend(list(reads(seq, dbg_cover=True, N=100)))
        random.shuffle(seq_reads)

        G = khmer.Nodegraph(K, 1e6, 4)
        sp = _oxli.StreamingPartitioner(G)

        for read in seq_reads:
            assert len(read) >= K
            sp.consume_sequence(read)
        assert sp.n_components == n_components

        comps = list(sp.components())
        comp = comps[0]
        assert len(comps) == n_components
        assert sp.n_components == (comp._n_created - comp._n_destroyed)
        assert sp.n_consumed == len(seq_reads)

    @pytest.mark.parametrize("n_components", list(range(1,101, 20)))
    @pytest.mark.parametrize("cov", [1,10,20])
    def test_write_components(self, random_sequence, cov, n_components, tmpdir):
        outfn = tmpdir.join('counts.csv')
        seqs = []
        for _ in range(n_components):
            seqs.append(random_sequence(exclude=''.join(seqs)))
        G = khmer.Countgraph(K, 1e6, 4)
        sp = _oxli.StreamingPartitioner(G)

        for seq in seqs:
            for _ in range(cov):
                sp.consume_sequence(seq)
        for seq in seqs:
            (med, _, _) = G.get_median_count(seq)
            assert med == cov
        assert sp.n_components == n_components

        sp.write_components(str(outfn))
        results = [line.strip().split(',') for line in outfn.open()]
        assert len(results) == n_components
        for row in results:
            assert abs(float(row[2])-float(cov)) < 2
