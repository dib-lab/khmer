import gc
import itertools
import random

from khmer import reverse_complement as revcomp
from khmer import reverse_hash as revhash
from khmer import forward_hash
from . import khmer_tst_utils as utils
from .graph_features import *

from khmer._oxli.graphlinks import Link, Junction, GraphLinker
from khmer import Nodegraph
import pytest


def teardown():
    utils.cleanup()


def test_get_junctions(right_tip_structure):
    graph, contig, L, HDN, R, tip = right_tip_structure
    linker = GraphLinker(graph)

    fw_choices, rc_choices = linker.get_junctions(contig)
    assert len(fw_choices) == 1
    assert len(rc_choices) == 0
    link = fw_choices.pop()

    assert link.u == forward_hash(HDN, K)
    assert link.v == forward_hash(R, K)


def test_get_junctions_rc(left_tip_structure):
    graph, contig, L, HDN, R, tip = left_tip_structure
    linker = GraphLinker(graph)

    fw_choices, rc_choices = linker.get_junctions(contig)
    assert len(fw_choices) == 0

    link = rc_choices.pop()

    assert link.u == forward_hash(HDN, K)
    assert link.v == forward_hash(L, K)




