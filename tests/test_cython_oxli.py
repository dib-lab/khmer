from __future__ import print_function
from __future__ import absolute_import

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

class TestStreamingPartitionerBasic:

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

    def test_merge_components(self, random_sequence):
        comp1 = random_sequence()
        comp2 = random_sequence(exclude=comp1)

        cg = khmer.Nodegraph(K, 1e5, 4)
        sp = _oxli.StreamingPartitioner(cg)
        
        sp.consume_sequence(comp1)
        sp.consume_sequence(comp2)
        assert sp.n_components == 2

        sp.consume_sequence(comp1 + comp2)
        assert sp.n_components == 1

        comps = list(sp.components())
        assert len(comps) == 1

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
