from __future__ import print_function
from __future__ import absolute_import

import gc
import itertools
import random

import khmer
from khmer._oxli.parsing import Sequence, FastxParser, BrokenPairedReader
from khmer._oxli.parsing import Alphabets, check_is_pair
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
        return fastx_fn
    return func


def test_FastxParser(create_fastx):
    reads = [Sequence('seq1/1', 'A' * 5),
              Sequence('seq1/2', 'A' * 4),
              Sequence('seq2/1', 'A' * 5),
              Sequence('seq3/1', 'A' * 3),
              Sequence('seq3/2', 'A' * 5)]
    x, m = gather(create_fastx(reads))

    expected = [('seq1/1', None),
                ('seq1/2', None),
                ('seq2/1', None),
                ('seq3/1', None),
                ('seq3/2', None)]

    assert x == expected, x

def test_FastxParser_sanitize(create_fastx):
    '''Test that A's are converted to N's when sanitize is True'''
    reads = [Sequence('seq1/1', 'N' * 5),
              Sequence('seq1/2', 'N' * 4)]
    parser = FastxParser(str(create_fastx(reads)), sanitize=True)

    parsed = [read for read in parser]
    assert parser.n_bad == 0
    assert parsed[0].sequence == 'A' * 5
    assert parsed[1].sequence == 'A' * 4

def test_FastxParser_no_sanitize(create_fastx):
    '''Test that N's remain when sanitize is False'''
    reads = [Sequence('seq1/1', 'N' * 5),
              Sequence('seq1/2', 'N' * 4)]
    parser = FastxParser(str(create_fastx(reads)), sanitize=False)

    parsed = [read for read in parser]
    assert parser.n_bad == 0
    assert parsed[0].sequence == 'N' * 5
    assert parsed[1].sequence == 'N' * 4

def test_FastxParser_on_invalid_sequence(create_fastx):
    '''Test that parser detects invalid sequence'''
    reads = [Sequence('seq1/1', 'XXX'),
              Sequence('seq1/2', 'A' * 4)]
    parser = FastxParser(str(create_fastx(reads)), sanitize=True)
    parsed = [read for read in parser]

    assert parser.n_bad == 1
    assert len(parsed) == 1
    assert parsed[0].sequence == 'A' * 4


def test_alphabet_wrapper():
    dna_simple = Alphabets.get('DNA_SIMPLE')
    assert len(dna_simple) == 4
    for b in 'ACGT':
        assert b in dna_simple

    with pytest.raises(ValueError):
        Alphabets.get('TEST')


def test_BrokenPairedReader_force_single(create_fastx):
    reads = [Sequence('seq1/1', 'A' * 5),
              Sequence('seq1/2', 'A' * 4),
              Sequence('seq2/1', 'A' * 5),
              Sequence('seq3/1', 'A' * 3),
              Sequence('seq3/2', 'A' * 5)]

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
    reads = [Sequence('seq1/1', 'A' * 5),
              Sequence('seq1/2', 'A' * 4),
              Sequence('seq3/1', 'A' * 3),
              Sequence('seq3/2', 'A' * 5)]


    x, n, m = gather_paired(create_fastx(reads), min_length=4, require_paired=True)

    expected = [('seq1/1', 'seq1/2')]
    assert x == expected, x
    assert m == 1
    assert n == 0, n


def test_BrokenPairedReader_OnPairs_filter_length_require_paired_2(create_fastx):
    reads = [Sequence('seq1/1', 'A' * 5),
              Sequence('seq1/2', 'A' * 4),
              Sequence('seq3/1', 'A' * 5),
              Sequence('seq3/2', 'A' * 3)]
    x, n, m = gather_paired(create_fastx(reads), min_length=4, require_paired=True)

    expected = [('seq1/1', 'seq1/2')]
    assert x == expected, x
    assert m == 1
    assert n == 0, n


def test_BrokenPairedReader_OnPairs_filter_length_required_paired_3(create_fastx):
    reads = [Sequence('seq1/1', 'A' * 5),
              Sequence('seq1/2', 'A' * 4),
              Sequence('seq3/1', 'A' * 3),
              Sequence('seq3/2', 'A' * 3)]
    x, n, m = gather_paired(create_fastx(reads), min_length=4, require_paired=True)

    expected = [('seq1/1', 'seq1/2')]
    assert x == expected, x
    assert m == 1
    assert n == 0, n


def test_BrokenPairedReader_OnPairs_filter_length_require_paired_4(create_fastx):
    reads = [Sequence('seq1/1', 'A' * 3),
              Sequence('seq1/2', 'A' * 3),
              Sequence('seq3/1', 'A' * 5),
              Sequence('seq3/2', 'A' * 5)]

    x, n, m = gather_paired(create_fastx(reads), min_length=4, require_paired=True)

    expected = [('seq3/1', 'seq3/2')]
    assert x == expected, x
    assert m == 1
    assert n == 0, n

'''
def test_BrokenPairedReader_lowercase():
    reads = [Sequence('seq1/1', 'acgtn'),
              Sequence('seq1/2', 'AcGtN'),
              Sequence('seq1/2', 'aCgTn')]

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


def test_check_is_pair_1():
    read1 = Sequence(name='seq', quality='###', sequence='AAA')
    read2 = Sequence(name='seq2', quality='###', sequence='AAA')

    assert not check_is_pair(read1, read2)


def test_check_is_pair_2():
    read1 = Sequence(name='seq/1', quality='###', sequence='AAA')
    read2 = Sequence(name='seq/2', quality='###', sequence='AAA')

    assert check_is_pair(read1, read2)


def test_check_is_pair_3_fq():
    read1 = Sequence(name='seq 1::', quality='###', sequence='AAA')
    read2 = Sequence(name='seq 2::', quality='###', sequence='AAA')

    assert check_is_pair(read1, read2)


def test_check_is_pair_3_broken_fq_1():
    read1 = Sequence(name='seq', quality='###', sequence='AAA')
    read2 = Sequence(name='seq 2::', quality='###', sequence='AAA')

    assert not check_is_pair(read1, read2)


def test_check_is_pair_3_broken_fq_2():
    read1 = Sequence(name='seq 1::', quality='###', sequence='AAA')
    read2 = Sequence(name='seq', quality='###', sequence='AAA')

    assert not check_is_pair(read1, read2)


def test_check_is_pair_3_fa():
    read1 = Sequence(name='seq 1::', sequence='AAA')
    read2 = Sequence(name='seq 2::', sequence='AAA')

    assert check_is_pair(read1, read2)


def test_check_is_pair_4():
    read1 = Sequence(name='seq/1', quality='###', sequence='AAA')
    read2 = Sequence(name='seq/2', sequence='AAA')

    try:
        check_is_pair(read1, read2)
        assert False                    # check_is_pair should fail here.
    except ValueError:
        pass


def test_check_is_pair_4b():
    read1 = Sequence(name='seq/1', sequence='AAA')
    read2 = Sequence(name='seq/2', quality='###', sequence='AAA')

    try:
        check_is_pair(read1, read2)
        assert False                    # check_is_pair should fail here.
    except ValueError:
        pass


def test_check_is_pair_5():
    read1 = Sequence(name='seq/1', sequence='AAA')
    read2 = Sequence(name='seq/2', sequence='AAA')

    assert check_is_pair(read1, read2)


def test_check_is_pair_6():
    read1 = Sequence(name='seq1', sequence='AAA')
    read2 = Sequence(name='seq2', sequence='AAA')

    assert not check_is_pair(read1, read2)


def test_check_is_pair_7():
    read1 = Sequence(name='seq/2', sequence='AAA')
    read2 = Sequence(name='seq/1', sequence='AAA')

    assert not check_is_pair(read1, read2)


def test_check_is_right():
    assert not check_is_right('seq1/1')
    assert not check_is_right('seq1 1::N')
    assert check_is_right('seq1/2')
    assert check_is_right('seq1 2::N')

    assert not check_is_right('seq')
    assert not check_is_right('seq 2')


def test_check_is_left():
    assert check_is_left('seq1/1')
    assert check_is_left('seq1 1::N')
    assert not check_is_left('seq1/2')
    assert not check_is_left('seq1 2::N')

    assert not check_is_left('seq')
    assert not check_is_left('seq 1')

    assert check_is_left(
        '@HWI-ST412:261:d15khacxx:8:1101:3149:2157 1:N:0:ATCACG')



