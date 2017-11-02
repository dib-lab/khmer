import gc
import itertools
import random

from khmer import reverse_complement as revcomp
from khmer import reverse_hash as revhash
from khmer import forward_hash
from . import khmer_tst_utils as utils
from .khmer_tst_utils import _equals_rc, _contains_rc
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
    assert len(nodes) == 1
    assert compactor.n_nodes == 1
    assert compactor.n_edges == 3

    node = nodes[0]
    assert _equals_rc(node.sequence, HDN)

    in_edges = list(node.in_edges())
    out_edges = list(node.out_edges())

    if len(in_edges) == 1:
        _, in_edge = in_edges[0]
        assert len(out_edges) == 2
        (_, edge_contig), (_, edge_tip) = out_edges 
        if len(edge_tip) > len(edge_contig):
            edge_contig, edge_tip = edge_tip, edge_contig
        assert _equals_rc(contig, in_edge.sequence + node.sequence +
                                  edge_contig.sequence)
    else:
        _, out_edge = out_edges[0]
        assert len(in_edges) == 2
        (_, edge_contig), (_, edge_tip) = in_edges
        if len(edge_tip) > len(edge_contig):
            edge_contig, edge_tip = edge_tip, edge_contig
        assert _equals_rc(contig, edge_contig.sequence + node.sequence +
                                  out_edge.sequence)

    for node in nodes:
        print(node)
        print('in edges:')
        for base, edge in node.in_edges():
            print(base, edge)
        
        print('out edges:')
        for base, edge in node.out_edges():
            print(base, edge)

    print("Contig FWD:", contig, len(contig))
    print("Contig RC:", revcomp(contig))
    print("HDN: ", repr(HDN))
    print("Tip FW:", tip, len(tip))
    print("Tip RC:", revcomp(tip))
    print("R FW:", R)
    print("R RC:", revcomp(R))

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
