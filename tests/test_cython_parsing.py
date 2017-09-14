
import gc
import itertools
import random

import khmer
from khmer._oxli.parsing import Sequence, FastxParser, SanitizedFastxParser
from khmer._oxli.parsing import BrokenPairedReader, Alphabets, check_is_pair
from khmer._oxli.parsing import check_is_right, check_is_left
from khmer.khmer_args import estimate_optimal_with_K_and_f as optimal_fp
from khmer import reverse_complement as revcomp
from khmer import reverse_hash as revhash
from . import khmer_tst_utils as utils

import pytest
import screed


def teardown():
    utils.cleanup()


@pytest.fixture
def create_fastx(tmpdir, as_str=True):
    def func(reads, fmt='fa'):
        assert fmt in ['fa', 'fq']
        fastx_fn = tmpdir.join('test.' + fmt)
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
        return str(fastx_fn) if as_str else fastx_fn
    return func


def sequences_eq(seqs_A, seqs_B):
    for seq1, seq2 in zip(seqs_A, seqs_B):
        if seq1 is None:
            assert seq1 is seq2
        else:
            assert seq1.name == seq2.name
            assert seq1.sequence == seq2.sequence


def test_FastxParser(create_fastx):
    expected = [Sequence('seq1/1', 'A' * 5),
                Sequence('seq1/2', 'A' * 4),
                Sequence('seq2/1', 'A' * 5),
                Sequence('seq3/1', 'A' * 3),
                Sequence('seq3/2', 'A' * 5)]
    parser = FastxParser(create_fastx(expected))
    result = list(parser)

    assert len(expected) == len(result)
    assert all((x == y) for x, y in zip(expected, result))


def test_SanitizedFastxParser_convert_Ns(create_fastx):
    '''Test that A's are converted to N's'''
    expected = [Sequence('seq1/1', 'N' * 5),
                Sequence('seq1/2', 'N' * 4)]
    parser = SanitizedFastxParser(create_fastx(expected),
                                  alphabet='DNAN_SIMPLE')
    result = list(parser)

    assert parser.n_bad == 0
    assert len(result) == 2
    assert result[0].sequence == 'A' * 5
    assert result[1].sequence == 'A' * 4


def test_SanitizedFastxParser_no_convert_Ns(create_fastx):
    expected = [Sequence('seq1/1', 'N' * 5),
                Sequence('seq1/2', 'N' * 4)]
    parser = SanitizedFastxParser(create_fastx(expected),
                                  alphabet='DNAN_SIMPLE',
                                  convert_n=False)
    result = list(parser)

    assert parser.n_bad == 0
    assert len(result) == 2
    assert result[0].sequence == 'N' * 5
    assert result[1].sequence == 'N' * 4


def test_SanitizedFastxParser_invalid(create_fastx):
    '''Test that parser detects invalid sequence'''
    expected = [Sequence('seq1/1', 'XXX'),
                Sequence('seq1/2', 'A' * 4)]
    parser = SanitizedFastxParser(create_fastx(expected))
    result = list(parser)

    assert parser.n_bad == 1
    assert len(result) == 1
    assert result[0].sequence == 'A' * 4


def test_SanitizedFastxParser_lowercase(create_fastx):
    reads = [Sequence('seq1/1', 'acgtn'),
             Sequence('seq1/2', 'AcGtN'),
             Sequence('seq1/2', 'aCgTn')]

    parser = SanitizedFastxParser(create_fastx(reads), convert_n=False)
    result = list(parser)

    assert result[0].sequence == 'ACGTN'
    assert result[1].sequence == 'ACGTN'
    assert result[2].sequence == 'ACGTN'


def test_alphabet_wrapper():
    dna_simple = Alphabets.get('DNA_SIMPLE')
    assert len(dna_simple) == 4
    for b in 'ACGT':
        assert b in dna_simple

    with pytest.raises(ValueError):
        Alphabets.get('TEST')


def gather_paired(stream, **kw):
    itr = BrokenPairedReader(stream, **kw)

    x = []
    m = 0
    num = 0
    for num, is_pair, read1, read2 in itr:
        x.append((read1.name if read1 is not None else None,
                  read2.name if read2 is not None else None))
        m += 1

    return x, num, m


class Test_BrokenPairedReader(object):
    reads = [Sequence(name='seq1/1', sequence='A' * 5),
             Sequence(name='seq1/2', sequence='A' * 4),
             Sequence(name='seq2/1', sequence='A' * 5),
             Sequence(name='seq3/1', sequence='A' * 3),
             Sequence(name='seq3/2', sequence='A' * 5)]

    @pytest.mark.parametrize("parser", [FastxParser, SanitizedFastxParser])
    def testDefault(self, parser, create_fastx):
        x, n, m = gather_paired(parser(create_fastx(self.reads)),
                                min_length=1)

        expected = [('seq1/1', 'seq1/2'),
                    ('seq2/1', None),
                    ('seq3/1', 'seq3/2')]
        assert x == expected, x
        assert m == 3
        assert n == 3, n

    @pytest.mark.parametrize("parser", [FastxParser, SanitizedFastxParser])
    def testMinLength(self, parser, create_fastx):
        x, n, m = gather_paired(parser(create_fastx(self.reads)),
                                min_length=3)

        expected = [('seq1/1', 'seq1/2'),
                    ('seq2/1', None),
                    ('seq3/1', 'seq3/2')]
        assert x == expected, x
        assert m == 3
        assert n == 3, n

    @pytest.mark.parametrize("parser", [FastxParser, SanitizedFastxParser])
    def testMinLength_2(self, parser, create_fastx):
        x, n, m = gather_paired(parser(create_fastx(self.reads)),
                                min_length=4)

        expected = [('seq1/1', 'seq1/2'),
                    ('seq2/1', None),
                    (None, 'seq3/2')]
        assert x == expected, x
        assert m == 3
        assert n == 3, n

    @pytest.mark.parametrize("parser", [FastxParser, SanitizedFastxParser])
    def testForceSingle(self, parser, create_fastx):
        x, n, m = gather_paired(parser(create_fastx(self.reads)),
                                force_single=True)

        expected = [('seq1/1', None),
                    ('seq1/2', None),
                    ('seq2/1', None),
                    ('seq3/1', None),
                    ('seq3/2', None)]
        assert x == expected, x
        assert m == 5
        assert n == 4, n

    @pytest.mark.parametrize("parser", [FastxParser, SanitizedFastxParser])
    def testForceSingleAndMinLength(self, parser, create_fastx):
        x, n, m = gather_paired(parser(create_fastx(self.reads)),
                                min_length=5, force_single=True)

        expected = [('seq1/1', None),
                    ('seq2/1', None),
                    ('seq3/2', None)]
        assert x == expected, x
        assert m == 3, m
        assert n == 2, n

    @pytest.mark.parametrize("parser", [FastxParser, SanitizedFastxParser])
    def testRequirePairedAndMinLength_HalfPass(self, parser, create_fastx):
        reads = [Sequence('seq1/1', 'A' * 5),
                 Sequence('seq1/2', 'A' * 4),
                 Sequence('seq3/1', 'A' * 3),
                 Sequence('seq3/2', 'A' * 5)]

        reader = BrokenPairedReader(parser(create_fastx(reads)),
                                    min_length=4, require_paired=True)

        result = []
        for n, paired, first, second in reader:
            result.append((first, second))

        assert len(result) == 1
        assert n == 0
        l, r = result[0]
        assert l == reads[0]
        assert r == reads[1]

    @pytest.mark.parametrize("parser", [FastxParser, SanitizedFastxParser])
    def testRequirePairedAndMinLength_SwappedHalfPass(self, parser,
                                                      create_fastx):
        reads = [Sequence('seq1/1', 'A' * 5),
                 Sequence('seq1/2', 'A' * 4),
                 Sequence('seq3/1', 'A' * 5),
                 Sequence('seq3/2', 'A' * 3)]

        reader = BrokenPairedReader(parser(create_fastx(reads)),
                                    min_length=4, require_paired=True)

        result = []
        for n, paired, first, second in reader:
            result.append((first, second))

        assert n == 0
        assert len(result) == 1
        l, r = result[0]
        assert l == reads[0]
        assert r == reads[1]

    @pytest.mark.parametrize("parser", [FastxParser, SanitizedFastxParser])
    def testRequirePairedAndMinLength_NeitherPass(self, parser, create_fastx):
        reads = [Sequence('seq1/1', 'A' * 5),
                 Sequence('seq1/2', 'A' * 4),
                 Sequence('seq3/1', 'A' * 3),
                 Sequence('seq3/2', 'A' * 3)]

        reader = BrokenPairedReader(parser(create_fastx(reads)),
                                    min_length=4, require_paired=True)

        result = []
        for n, paired, first, second in reader:
            result.append((first, second))

        assert n == 0
        assert len(result) == 1
        l, r = result[0]
        assert l == reads[0]
        assert r == reads[1]

    @pytest.mark.parametrize("parser", [FastxParser, SanitizedFastxParser])
    def testRequirePairedAndMinLength_SwappedNeitherPass(self, parser,
                                                         create_fastx):
        reads = [Sequence('seq1/1', 'A' * 3),
                 Sequence('seq1/2', 'A' * 3),
                 Sequence('seq3/1', 'A' * 5),
                 Sequence('seq3/2', 'A' * 5)]

        reader = BrokenPairedReader(parser(create_fastx(reads)),
                                    min_length=4, require_paired=True)

        result = []
        for n, paired, first, second in reader:
            result.append((first, second))

        assert n == 0
        assert len(result) == 1
        l, r = result[0]
        assert l == reads[2]
        assert r == reads[3]


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


class Test_Sequence(object):

    name = 'Test'
    sequence = 'ACGT'
    quality = '####'
    description = 'The nucleotides'
    cleaned = 'aaaa'

    def test_init_name_and_sequence(self):
        s = Sequence(name=self.name, sequence=self.sequence)
        assert s.name == self.name
        assert s.sequence == self.sequence
        assert s.quality is None
        assert s.description is None
        assert s.cleaned_seq == self.sequence

    def test_init_name_only(self):
        s = Sequence(name=self.name)
        assert s.name is None
        assert s.sequence is None
        assert s.quality is None
        assert s.description is None
        assert s.cleaned_seq is None

    def test_init_sequence_only(self):
        s = Sequence(sequence=self.sequence)
        assert s.name is None
        assert s.sequence is None
        assert s.quality is None
        assert s.description is None
        assert s.cleaned_seq is None

    def test_init_with_cleaned_seq(self):
        s = Sequence(name=self.name, sequence=self.sequence,
                     cleaned_seq=self.cleaned)
        assert s.name == self.name
        assert s.sequence == self.sequence
        assert s.quality is None
        assert s.description is None
        assert s.cleaned_seq == self.cleaned
