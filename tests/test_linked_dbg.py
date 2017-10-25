import gc
import itertools
import random

from khmer import reverse_complement as revcomp
from khmer import reverse_hash as revhash
from khmer import forward_hash
from . import khmer_tst_utils as utils
from .graph_features import *

from khmer._oxli.graphlinks import StreamingCompactor
from khmer import Nodegraph
import pytest


def teardown():
    utils.cleanup()


def test_compact_fork(right_tip_structure):
    '''Should have no links. Need two junctions.
    '''
    graph, contig, L, HDN, R, tip = right_tip_structure
    compactor = StreamingCompactor(graph)
    compactor.update(contig)
    compactor.report()

    nodes = list(compactor.sequence_nodes(contig))
    assert len(nodes) == 0
    assert compactor.n_nodes == 0
    assert compactor.n_edges == 0



    #assert junction['u'] == forward_hash(L, K)
    #assert junction['v'] == forward_hash(HDN, K)
    #assert junction['w'] == forward_hash(R, K)

'''
def test_links_bubble(snp_bubble_structure):
    graph, wildtype_sequence, _, HDN_L, HDN_R = snp_bubble_structure
    linker = GraphLinker(graph)

    # thread the wildtype half of the bubble
    linker.add_links(wildtype_sequence)
    linker.report()


    links = list(linker.get_links(wildtype_sequence))
    assert len(links) == 2
    
    link_a, link_b = links
    if link_a[0]['v'] == forward_hash(HDN_L, K):
        assert link_b[0]['v'] == forward_hash(HDN_R, K)
    elif link_a[0]['v'] == forward_hash(HDN_R, K):
        assert link_b[0]['v'] == forward_hash(HDN_L, K)
'''
