
import screed
import khmer
import os
import sys
import pytest
from . import khmer_tst_utils as utils
from .graph_features import random_sequence, kmers, K

from khmer._oxli.hashing import SequenceHasher
from khmer import forward_hash


def test_sequence_hasher(random_sequence):
    sequence = random_sequence()

    hashvals = []
    for h in SequenceHasher(sequence, K):
        hashvals.append(h)
    exp_hashvals = [forward_hash(kmer, K) for kmer in kmers(sequence)]

    assert len(hashvals) == len(exp_hashvals)
    assert hashvals == exp_hashvals
