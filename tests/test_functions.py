#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import khmer
from nose.plugins.attrib import attr
import os
import khmer_tst_utils as utils
import collections
from khmer.utils import check_is_pair


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
        print ksize, table_size, n_tables

        assert(ksize) == 25
        assert table_size == size
        assert n_tables == 4

        try:
            os.remove(fn)
        except OSError as e:
            print >>sys.stder, '...failed to remove {fn}'.format(fn)


def test_extract_hashbits_info():
    fn = utils.get_temp_filename('test_extract_hashbits.pt')
    for size in [1e6, 2e6, 5e6, 1e7]:
        ht = khmer.Hashbits(25, size, 4)
        ht.save(fn)

        info = khmer.extract_hashbits_info(fn)
        ksize, table_size, n_tables, _, _ = info
        print ksize, table_size, n_tables

        assert(ksize) == 25
        assert table_size == size
        assert n_tables == 4

        try:
            os.remove(fn)
        except OSError as e:
            print >>sys.stderr, '...failed to remove {fn}'.format(fn)


FakeFQRead = collections.namedtuple('Read', ['name', 'accuracy', 'sequence'])
FakeFastaRead = collections.namedtuple('Read', ['name', 'sequence'])


def test_check_is_pair_1():
    read1 = FakeFQRead(name='seq', accuracy='###', sequence='AAA')
    read2 = FakeFQRead(name='seq2', accuracy='###', sequence='AAA')

    assert not check_is_pair(read1, read2)


def test_check_is_pair_2():
    read1 = FakeFQRead(name='seq/1', accuracy='###', sequence='AAA')
    read2 = FakeFQRead(name='seq/2', accuracy='###', sequence='AAA')

    assert check_is_pair(read1, read2)


def test_check_is_pair_3():
    read1 = FakeFQRead(name='seq 1::', accuracy='###', sequence='AAA')
    read2 = FakeFQRead(name='seq 2::', accuracy='###', sequence='AAA')

    assert check_is_pair(read1, read2)


def test_check_is_pair_4():
    read1 = FakeFQRead(name='seq/1', accuracy='###', sequence='AAA')
    read2 = FakeFastaRead(name='seq/2', sequence='AAA')

    try:
        check_is_pair(read1, read2)
        assert False                    # check_is_pair should fail here.
    except ValueError:
        pass


def test_check_is_pair_4b():
    read1 = FakeFastaRead(name='seq/1', sequence='AAA')
    read2 = FakeFQRead(name='seq/2', accuracy='###', sequence='AAA')

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
