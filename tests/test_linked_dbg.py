import gc
import itertools
import random

from khmer import reverse_complement as revcomp
from khmer import reverse_hash as revhash
from khmer import forward_hash
from . import khmer_tst_utils as utils
from .graph_features import *

from khmer._oxli.graphlinks import Link, LinkPath, GraphLinker
from khmer import Nodegraph
import pytest


def teardown():
    utils.cleanup()


def test_get_junction_choices(left_tip_structure):
    graph, contig, L, HDN, R, tip = left_tip_structure
    linker = GraphLinker(graph)

    choices = linker.get_junction_choices(contig)
    assert len(choices) == 1
    link = choices.pop()

    assert link.u == forward_hash(L, K)
    assert link.v == forward_hash(HDN, K)


def test_get_junction_choices_rc(left_tip_structure):
    graph, contig, L, HDN, R, tip = left_tip_structure
    linker = GraphLinker(graph)

    choices = linker.get_junction_choices(revcomp(contig))
    assert len(choices) == 0


def test_get_junction_choices_right_rc(right_tip_structure):
    graph, contig, L, HDN, R, tip = right_tip_structure
    linker = GraphLinker(graph)

    choices = linker.get_junction_choices(revcomp(contig))
    assert len(choices) == 1
    link = choices.pop()

    assert link.u == forward_hash(R, K)
    assert link.v == forward_hash(HDN, K)


