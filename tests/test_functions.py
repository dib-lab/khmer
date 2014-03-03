#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import khmer
from nose.plugins.attrib import attr
import os


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


def test_get_primes():
    primes = khmer.get_n_primes_near_x(7, 20)

    assert primes == [19, 17, 13, 11, 7, 5, 3]


def test_extract_countinghash_info():
    fn = 'test_extract_counting.ht'
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
    fn = 'test_extract_hashbits.ht'
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
