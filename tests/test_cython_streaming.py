from __future__ import print_function
from __future__ import absolute_import

import gc
import itertools
import random

import khmer
from khmer._oxli.streaming import Sequence, FastxParser, BrokenPairedReader
from khmer.khmer_args import estimate_optimal_with_K_and_f as optimal_fp
from khmer import reverse_complement as revcomp
from khmer import reverse_hash as revhash
from . import khmer_tst_utils as utils

import pytest
import screed


def teardown():
    utils.cleanup()

def gather(filename, **kw):
    stream = FastxParser(str(filename))

    x = []
    m = 0
    for read in stream:
        x.append((read.name, None))
        m += 1

    return x,m

def gather_paired(filename, **kw):
    stream = FastxParser(str(filename))
    itr = BrokenPairedReader(stream, **kw)

    x = []
    m = 0
    num = 0
    for num, is_pair, read1, read2 in itr:
        if is_pair:
            x.append((read1.name, read2.name))
        else:
            x.append((read1.name, None))
        m += 1

    return x, num, m


@pytest.fixture
def create_fastx(tmpdir):
    def func(reads, fmt='fa'):
        assert fmt in ['fa','fq']
        fastx_fn = tmpdir.join('test.'+fmt)
        for record in reads:
            if fmt == 'fa':
                fastx_fn.write('>{0}\n{1}\n'.format(record.name,
                                                    record.sequence),
                               mode='a')
            else:
                fastx_fn.write('@{0}\n{1}\n+\n{2}\n'.format(record.name,
                                                        record.sequence,
                                                        record.quality),
                               mode='a')
        print(fastx_fn.read())
        return fastx_fn
    return func


def test_FastxParser(create_fastx):
    reads = [Sequence.new('seq1/1', 'A' * 5),
              Sequence.new('seq1/2', 'A' * 4),
              Sequence.new('seq2/1', 'A' * 5),
              Sequence.new('seq3/1', 'A' * 3),
              Sequence.new('seq3/2', 'A' * 5)]
    x, m = gather(create_fastx(reads))

    expected = [('seq1/1', None),
                ('seq1/2', None),
                ('seq2/1', None),
                ('seq3/1', None),
                ('seq3/2', None)]

    assert x == expected, x

def test_BrokenPairedReader_force_single(create_fastx):
    reads = [Sequence.new('seq1/1', 'A' * 5),
              Sequence.new('seq1/2', 'A' * 4),
              Sequence.new('seq2/1', 'A' * 5),
              Sequence.new('seq3/1', 'A' * 3),
              Sequence.new('seq3/2', 'A' * 5)]

    x, n, m = gather_paired(create_fastx(reads), force_single=True)

    expected = [('seq1/1', None),
                ('seq1/2', None),
                ('seq2/1', None),
                ('seq3/1', None),
                ('seq3/2', None)]
    assert x == expected, x
    assert m == 5
    assert n == 4, n


def test_BrokenPairedReader_OnPairs_filter_length_require_paired(create_fastx):
    reads = [Sequence.new('seq1/1', 'A' * 5),
              Sequence.new('seq1/2', 'A' * 4),
              Sequence.new('seq3/1', 'A' * 3),
              Sequence.new('seq3/2', 'A' * 5)]


    x, n, m = gather_paired(create_fastx(reads), min_length=4, require_paired=True)

    expected = [('seq1/1', 'seq1/2')]
    assert x == expected, x
    assert m == 1
    assert n == 0, n


def test_BrokenPairedReader_OnPairs_filter_length_require_paired_2(create_fastx):
    reads = [Sequence.new('seq1/1', 'A' * 5),
              Sequence.new('seq1/2', 'A' * 4),
              Sequence.new('seq3/1', 'A' * 5),
              Sequence.new('seq3/2', 'A' * 3)]
    x, n, m = gather_paired(create_fastx(reads), min_length=4, require_paired=True)

    expected = [('seq1/1', 'seq1/2')]
    assert x == expected, x
    assert m == 1
    assert n == 0, n


def test_BrokenPairedReader_OnPairs_filter_length_required_paired_3(create_fastx):
    reads = [Sequence.new('seq1/1', 'A' * 5),
              Sequence.new('seq1/2', 'A' * 4),
              Sequence.new('seq3/1', 'A' * 3),
              Sequence.new('seq3/2', 'A' * 3)]
    x, n, m = gather_paired(create_fastx(reads), min_length=4, require_paired=True)

    expected = [('seq1/1', 'seq1/2')]
    assert x == expected, x
    assert m == 1
    assert n == 0, n


def test_BrokenPairedReader_OnPairs_filter_length_require_paired_4(create_fastx):
    reads = [Sequence.new('seq1/1', 'A' * 3),
              Sequence.new('seq1/2', 'A' * 3),
              Sequence.new('seq3/1', 'A' * 5),
              Sequence.new('seq3/2', 'A' * 5)]

    x, n, m = gather_paired(create_fastx(reads), min_length=4, require_paired=True)

    expected = [('seq3/1', 'seq3/2')]
    assert x == expected, x
    assert m == 1
    assert n == 0, n

'''
def test_BrokenPairedReader_lowercase():
    reads = [Sequence.new('seq1/1', 'acgtn'),
              Sequence.new('seq1/2', 'AcGtN'),
              Sequence.new('seq1/2', 'aCgTn')]

    results = []
    parser = FastxParser(create_fastx(reads))
    for num, is_pair, read1, read2 in broken_paired_reader(parser):
        results.append((read1, read2))

    a, b = results[0]
    assert a.sequence == 'acgtn'
    assert a.cleaned_seq == 'ACGTA'
    assert b.sequence == 'AcGtN'
    assert b.cleaned_seq == 'ACGTA'

    c, d = results[1]
    assert c.sequence == 'aCgTn'
    assert c.cleaned_seq == 'ACGTA'
    assert d is None
'''
