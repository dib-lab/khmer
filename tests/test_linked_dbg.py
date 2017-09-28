import gc
import itertools
import random

from khmer import reverse_complement as revcomp
from khmer import reverse_hash as revhash
from khmer import forward_hash
from . import khmer_tst_utils as utils
from .graph_features import *

from khmer._oxli.graphlinks import GraphLinker
from khmer import Nodegraph
import pytest


def teardown():
    utils.cleanup()


def test_add_links_one_junction(right_tip_structure):
    '''Should have no links. Need two junctions.
    '''
    graph, contig, L, HDN, R, tip = right_tip_structure
    linker = GraphLinker(graph)

    linker.add_links(contig)
    linker.report()

    links = list(linker.get_links(contig))
    assert len(links) == 0
    
    junctions = list(linker.get_junctions(contig))
    assert len(junctions) == 1

    print(junctions)



