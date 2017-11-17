import gc
import itertools
import random

from khmer import reverse_complement as revcomp
from khmer import reverse_hash as revhash
from khmer import forward_hash
from . import khmer_tst_utils as utils
from .khmer_tst_utils import _equals_rc, _contains_rc
from .graph_structure_fixtures import *

from khmer._oxli.graphlinks import StreamingCompactor
from khmer import Nodegraph
import pytest


def teardown():
    utils.cleanup()


def compare_tip_with_cdbg(rts, compactor):
    graph, contig, L, HDN, R, tip = rts

    nodes = list(compactor.sequence_nodes(contig))
    assert len(nodes) == 1

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
        #assert _equals_rc(contig, in_edge.sequence[:-K+1] + node.sequence +
        #                  edge_contig.sequence[K-1:])
    else:
        _, out_edge = out_edges[0]
        assert len(in_edges) == 2
        (_, edge_contig), (_, edge_tip) = in_edges
        if len(edge_tip) > len(edge_contig):
            edge_contig, edge_tip = edge_tip, edge_contig
        #assert _equals_rc(contig, edge_contig.sequence[:-K+1] + node.sequence +
        #                  out_edge.sequence[K-1:])


@using_ksize([21,25,31])
def test_compact_tip(ksize, right_tip_structure):
    '''Should have no links. Need two junctions.
    '''
    right_tip_structure = right_tip_structure()
    graph, contig, L, HDN, R, tip = right_tip_structure

    compactor = StreamingCompactor(graph)
    print(compactor.update(contig), 'cDBG updates...')
    compactor.report()

    compare_tip_with_cdbg(right_tip_structure, compactor)

    assert compactor.n_nodes == 1
    assert compactor.n_edges == 3

    for node in compactor.sequence_nodes(contig):
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


def test_compact_tip_double_update(right_tip_structure):
    right_tip_structure = right_tip_structure()
    graph, contig, L, HDN, R, tip = right_tip_structure

    compactor = StreamingCompactor(graph)
    print(compactor.update(contig), 'cDBG updates...')
    compactor.report()
    print(compactor.update(contig), 'cDBG updates...')
    compactor.report()

    compare_tip_with_cdbg(right_tip_structure, compactor)
    assert compactor.n_nodes == 1
    assert compactor.n_edges == 3


def test_compact_tip_revcomp_update(right_tip_structure):
    right_tip_structure = right_tip_structure()
    graph, contig, L, HDN, R, tip = right_tip_structure

    compactor = StreamingCompactor(graph)
    print(compactor.update(contig), 'cDBG updates...')
    compactor.report()

    print(compactor.update(revcomp(contig)), 'cDBG updates...')
    compactor.report()

    compare_tip_with_cdbg(right_tip_structure, compactor)
    assert compactor.n_nodes == 1
    assert compactor.n_edges == 3


def test_compact_two_tip_islands(left_tip_structure, right_tip_structure):
    right_tip_structure = right_tip_structure()
    graph, contig_r, L_r, HDN_r, R_r, tip_r = right_tip_structure
    left_tip_structure = left_tip_structure()
    _, contig_l, L_l, HDN_l, R_l, tip_l = left_tip_structure
    
    compactor = StreamingCompactor(graph)
    print(compactor.update(contig_l), 'cDBG updates from left')
    compactor.report()
    compare_tip_with_cdbg(left_tip_structure, compactor)
    assert compactor.n_nodes == 1
    assert compactor.n_edges == 3

    print(compactor.update(contig_r), 'cDBG updates from right')
    compactor.report()
    compare_tip_with_cdbg(right_tip_structure, compactor)
    assert compactor.n_nodes == 2
    assert compactor.n_edges == 6


def test_compact_tip_x_merge(left_tip_structure, right_tip_structure):
    right_tip_structure = right_tip_structure()
    graph, contig_r, L_r, HDN_r, R_r, tip_r = right_tip_structure
    left_tip_structure = left_tip_structure()
    _, contig_l, L_l, HDN_l, R_l, tip_l = left_tip_structure
    
    contig_merge = contig_l + contig_r
    graph.reset()
    
    compactor = StreamingCompactor(graph)
    compactor.consume(str(tip_l))
    print(compactor.consume_and_update(contig_l),
          'cDBG updates from left')
    compactor.report()
    compare_tip_with_cdbg(left_tip_structure, compactor)
    assert compactor.n_nodes == 1
    assert compactor.n_edges == 3

    compactor.consume(str(tip_r))
    print(compactor.consume_and_update(contig_merge), 
          'cDBG updates from right merge')
    compactor.report()
    compare_tip_with_cdbg(right_tip_structure, compactor)
    assert compactor.n_nodes == 2
    assert compactor.n_edges == 5


@using_ksize([21, 31])
def test_compact_triple_fork(right_triple_fork_structure):
    right_triple_fork_structure = right_triple_fork_structure()
    graph, core, L, HDN, R, top, bottom = right_triple_fork_structure

    compactor = StreamingCompactor(graph)
    compactor.update(core)
    compactor.report()

    assert compactor.n_nodes == 1
    assert compactor.n_edges == 4


@pytest.mark.parametrize('random_sequence', [100, 200], indirect=True)
def test_compact_trivial_edge(tandem_triple_forks, ksize):
    ttf = tandem_triple_forks()
    graph, core, L, HDN_l, HDN_r, R, top_l, bottom_l, top_r, bottom_r = ttf

    print('Core:', core[HDN_l.pos:], '\nHDN_l:', HDN_l, '\nHDN_r:', HDN_r,
          '\ntop_l:', top_l[HDN_l.pos:HDN_l.pos+2*ksize], 
          '\nbottom_l:', bottom_l[HDN_l.pos:HDN_l.pos+2*ksize],
          '\ntop_r:', top_r[HDN_r.pos:HDN_r.pos+2*ksize],
          '\nbottom_r:', bottom_r[HDN_r.pos:HDN_r.pos+2*ksize])
    br = '=' * 20
    graph.reset()
    compactor = StreamingCompactor(graph)
    print(br, 'ADD CORE', br) 
    compactor.consume_and_update(core)
    assert compactor.n_nodes == 0

    print(br, 'ADD top_l', br)
    compactor.consume_and_update(top_l)
    assert compactor.n_nodes == 1
    assert compactor.n_edges == 3

    print(br, 'ADD bottom_l', br)
    compactor.consume_and_update(bottom_l)
    assert compactor.n_nodes == 1
    assert compactor.n_edges == 4

    print(br, 'ADD top_r', br)
    compactor.consume_and_update(top_r)
    assert compactor.n_nodes == 2
    assert compactor.n_edges == 6

    print(br, 'ADD bottom_r', br)
    compactor.consume_and_update(bottom_r)
    assert compactor.n_nodes == 2
    assert compactor.n_edges == 7

    nodes = list(compactor.sequence_nodes(core))
    node_1, node_2 = nodes
    trivial, node_2_out = list(node_2.in_edges()), list(node_2.out_edges())
    if len(trivial) != 1:
        trivial, node_2_out = node_2_out, trivial
    _, trivial = trivial[0]

    assert trivial.edge_type == 'TRIVIAL'
    assert len(trivial) == ksize + 1

    assert HDN_l in trivial.sequence
    assert HDN_r in trivial.sequence
    assert node_1.degree == 4
    assert node_2.degree == 4


def test_compact_tip_split_merge():
    pass
