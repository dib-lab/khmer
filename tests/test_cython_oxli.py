from __future__ import print_function
from __future__ import absolute_import

import itertools
import random

import khmer
from khmer import _oxli
from khmer.khmer_args import estimate_optimal_with_K_and_f as optimal_fp
from khmer import reverse_complement as revcomp
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

    def test_get_nearest_component(self, random_sequence):
        comp1 =  random_sequence()
        comp2 = random_sequence(exclude=comp1)

        cg = khmer.Nodegraph(K, 1e5, 4)
        sp = _oxli.StreamingPartitioner(cg)
        assert False
        sp.consume_sequence(comp1)
        sp.consume_sequence(comp2)

        c = sp.get_nearest_component(comp1[:K])
        assert c.component_id == 0
