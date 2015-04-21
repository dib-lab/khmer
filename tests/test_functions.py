from __future__ import print_function
from __future__ import absolute_import
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
import khmer
from nose.plugins.attrib import attr
import os
from . import khmer_tst_utils as utils
import collections
from khmer.utils import (check_is_pair, broken_paired_reader, check_is_left,
                         check_is_right)
from khmer.kfile import check_input_files


def test_forward_hash():
    assert khmer.forward_hash('AAAA', 4) == 0
    assert khmer.forward_hash('TTTT', 4) == 0
    assert khmer.forward_hash('CCCC', 4) == 170
    assert khmer.forward_hash('GGGG', 4) == 170


def test_forward_hash_no_rc():
    h = khmer.forward_hash_no_rc('AAAA', 4)
    assert h == 0, h

    h = khmer.forward_hash_no_rc('TTTT', 4)
    assert h == 85, h

    h = khmer.forward_hash_no_rc('CCCC', 4)
    assert h == 170, h

    h = khmer.forward_hash_no_rc('GGGG', 4)
    assert h == 255, h


def test_reverse_hash():
    s = khmer.reverse_hash(0, 4)
    assert s == "AAAA"

    s = khmer.reverse_hash(85, 4)
    assert s == "TTTT"

    s = khmer.reverse_hash(170, 4)
    assert s == "CCCC"

    s = khmer.reverse_hash(255, 4)
    assert s == "GGGG"


def test_hash_murmur3():
    assert khmer.hash_murmur3('AAAA') == 526240128537019279
    assert khmer.hash_murmur3('TTTT') == 526240128537019279
    assert khmer.hash_murmur3('CCCC') == 14391997331386449225
    assert khmer.hash_murmur3('GGGG') == 14391997331386449225


def test_hash_no_rc_murmur3():
    h = khmer.hash_no_rc_murmur3('AAAA')
    assert h == 5231866503566620412, h

    h = khmer.hash_no_rc_murmur3('TTTT')
    assert h == 5753003579327329651, h

    h = khmer.hash_no_rc_murmur3('CCCC')
    assert h == 3789793362494378039, h

    h = khmer.hash_no_rc_murmur3('GGGG')
    assert h == 17519752047064575358, h


def test_get_primes():
    primes = khmer.get_n_primes_near_x(7, 20)

    assert primes == [19, 17, 13, 11, 7, 5, 3]


def test_extract_countinghash_info():
    fn = utils.get_temp_filename('test_extract_counting.ct')
    for size in [1e6, 2e6, 5e6, 1e7]:
        ht = khmer.new_counting_hash(25, size, 4)
        ht.save(fn)

        info = khmer.extract_countinghash_info(fn)
        ksize, table_size, n_tables, _, _, _ = info
        print(ksize, table_size, n_tables)

        assert(ksize) == 25
        assert table_size == size
        assert n_tables == 4

        try:
            os.remove(fn)
        except OSError as e:
            print('...failed to remove {fn}'.format(fn), file=sys.stder)


def test_extract_hashbits_info():
    fn = utils.get_temp_filename('test_extract_hashbits.pt')
    for size in [1e6, 2e6, 5e6, 1e7]:
        ht = khmer.Hashbits(25, size, 4)
        ht.save(fn)

        info = khmer.extract_hashbits_info(fn)
        ksize, table_size, n_tables, _, _ = info
        print(ksize, table_size, n_tables)

        assert(ksize) == 25
        assert table_size == size
        assert n_tables == 4

        try:
            os.remove(fn)
        except OSError as e:
            print('...failed to remove {fn}'.format(fn), file=sys.stderr)


def test_check_file_status_kfile():
    fn = utils.get_temp_filename('thisfiledoesnotexist')
    check_file_status_exited = False
    try:
        check_input_files(fn, False)
    except SystemExit:
        check_file_status_exited = True
    assert check_file_status_exited


def test_check_file_status_kfile_force():
    fn = utils.get_temp_filename('thisfiledoesnotexist')
    try:
        check_input_files(fn, True)
    except OSError as e:
        assert False


FakeFQRead = collections.namedtuple('Read', ['name', 'quality', 'sequence'])
FakeFastaRead = collections.namedtuple('Read', ['name', 'sequence'])


def test_check_is_pair_1():
    read1 = FakeFQRead(name='seq', quality='###', sequence='AAA')
    read2 = FakeFQRead(name='seq2', quality='###', sequence='AAA')

    assert not check_is_pair(read1, read2)


def test_check_is_pair_2():
    read1 = FakeFQRead(name='seq/1', quality='###', sequence='AAA')
    read2 = FakeFQRead(name='seq/2', quality='###', sequence='AAA')

    assert check_is_pair(read1, read2)


def test_check_is_pair_3_fq():
    read1 = FakeFQRead(name='seq 1::', quality='###', sequence='AAA')
    read2 = FakeFQRead(name='seq 2::', quality='###', sequence='AAA')

    assert check_is_pair(read1, read2)


def test_check_is_pair_3_broken_fq_1():
    read1 = FakeFQRead(name='seq', quality='###', sequence='AAA')
    read2 = FakeFQRead(name='seq 2::', quality='###', sequence='AAA')

    assert not check_is_pair(read1, read2)


def test_check_is_pair_3_broken_fq_2():
    read1 = FakeFQRead(name='seq 1::', quality='###', sequence='AAA')
    read2 = FakeFQRead(name='seq', quality='###', sequence='AAA')

    assert not check_is_pair(read1, read2)


def test_check_is_pair_3_fa():
    read1 = FakeFastaRead(name='seq 1::', sequence='AAA')
    read2 = FakeFastaRead(name='seq 2::', sequence='AAA')

    assert check_is_pair(read1, read2)


def test_check_is_pair_4():
    read1 = FakeFQRead(name='seq/1', quality='###', sequence='AAA')
    read2 = FakeFastaRead(name='seq/2', sequence='AAA')

    try:
        check_is_pair(read1, read2)
        assert False                    # check_is_pair should fail here.
    except ValueError:
        pass


def test_check_is_pair_4b():
    read1 = FakeFastaRead(name='seq/1', sequence='AAA')
    read2 = FakeFQRead(name='seq/2', quality='###', sequence='AAA')

    try:
        check_is_pair(read1, read2)
        assert False                    # check_is_pair should fail here.
    except ValueError:
        pass


def test_check_is_pair_5():
    read1 = FakeFastaRead(name='seq/1', sequence='AAA')
    read2 = FakeFastaRead(name='seq/2', sequence='AAA')

    assert check_is_pair(read1, read2)


def test_check_is_pair_6():
    read1 = FakeFastaRead(name='seq1', sequence='AAA')
    read2 = FakeFastaRead(name='seq2', sequence='AAA')

    assert not check_is_pair(read1, read2)


def test_check_is_pair_7():
    read1 = FakeFastaRead(name='seq/2', sequence='AAA')
    read2 = FakeFastaRead(name='seq/1', sequence='AAA')

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


class Test_BrokenPairedReader(object):
    stream = [FakeFastaRead(name='seq1/1', sequence='A' * 5),
              FakeFastaRead(name='seq1/2', sequence='A' * 4),
              FakeFastaRead(name='seq2/1', sequence='A' * 5),
              FakeFastaRead(name='seq3/1', sequence='A' * 3),
              FakeFastaRead(name='seq3/2', sequence='A' * 5)]

    def gather(self, **kw):
        iter = broken_paired_reader(self.stream, **kw)

        x = []
        m = 0
        for n, is_pair, read1, read2 in iter:
            if is_pair:
                x.append((read1.name, read2.name))
            else:
                x.append((read1.name, None))
            m += 1

        return x, n, m

    def testDefault(self):
        x, n, m = self.gather(min_length=1)

        expected = [('seq1/1', 'seq1/2'),
                    ('seq2/1', None),
                    ('seq3/1', 'seq3/2')]
        assert x == expected, x
        assert m == 3
        assert n == 3, n

    def testMinLength(self):
        x, n, m = self.gather(min_length=3)

        expected = [('seq1/1', 'seq1/2'),
                    ('seq2/1', None),
                    ('seq3/1', 'seq3/2')]
        assert x == expected, x
        assert m == 3
        assert n == 3, n

    def testMinLength_2(self):
        x, n, m = self.gather(min_length=4)

        expected = [('seq1/1', 'seq1/2'),
                    ('seq2/1', None),
                    ('seq3/2', None)]
        assert x == expected, x
        assert m == 3
        assert n == 3, n

    def testForceSingle(self):
        x, n, m = self.gather(force_single=True)

        expected = [('seq1/1', None),
                    ('seq1/2', None),
                    ('seq2/1', None),
                    ('seq3/1', None),
                    ('seq3/2', None)]
        assert x == expected, x
        assert m == 5
        assert n == 4, n

    def testForceSingleAndMinLength(self):
        x, n, m = self.gather(min_length=5, force_single=True)

        expected = [('seq1/1', None),
                    ('seq2/1', None),
                    ('seq3/2', None)]
        assert x == expected, x
        assert m == 3, m
        assert n == 2, n
