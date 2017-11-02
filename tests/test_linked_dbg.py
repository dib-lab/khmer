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


def compare_right_tip_with_cdbg(rts, compactor):
    graph, contig, L, HDN, R, tip = rts

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


def test_compact_fork(right_tip_structure):
    '''Should have no links. Need two junctions.
    '''
    graph, contig, L, HDN, R, tip = right_tip_structure

    compactor = StreamingCompactor(graph)
    print(compactor.update(contig), 'cDBG updates...')
    compactor.report()

    compare_right_tip_with_cdbg(right_tip_structure, compactor)

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


def test_compact_fork_double_update(right_tip_structure):
    graph, contig, L, HDN, R, tip = right_tip_structure

    compactor = StreamingCompactor(graph)
    print(compactor.update(contig), 'cDBG updates...')
    compactor.report()
    print(compactor.update(contig), 'cDBG updates...')
    compactor.report()

    compare_right_tip_with_cdbg(right_tip_structure, compactor)


def test_compact_fork_revcomp_update(right_tip_structure):
    graph, contig, L, HDN, R, tip = right_tip_structure

    compactor = StreamingCompactor(graph)
    print(compactor.update(contig), 'cDBG updates...')
    compactor.report()

    print(compactor.update(revcomp(contig)), 'cDBG updates...')
    compactor.report()

    compare_right_tip_with_cdbg(right_tip_structure, compactor)

