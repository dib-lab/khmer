from __future__ import print_function
from __future__ import absolute_import, unicode_literals
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# pylint: disable=missing-docstring,protected-access
import gzip

import os
import shutil

import khmer
from . import khmer_tst_utils as utils
from khmer import ReadParser
import screed

import nose
from nose.plugins.attrib import attr

MAX_COUNT = 255
MAX_BIGCOUNT = 65535

#

# from http://www.rsok.com/~jrm/printprimes.html
PRIMES_1m = [1000003, 1009837]
PRIMES_100m = [100009979, 100000007]
PRIMES_1b = [1000000007, 1000000919]
PRIMES_2b = [1999999973, 1999999943]
PRIMES_4b = [4000000007, 4000000009]
PRIMES_8b = [8000000011, 8000000051]

DNA = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"


def teardown():
    utils.cleanup()


class Test_CountingHash(object):

    def setup(self):
        self.hi = khmer.CountingHash(12, PRIMES_1m)

    def test_failed_get(self):
        GG = 'G' * 12                   # forward_hash: 11184810
        GGhash = khmer.forward_hash(GG, 12)
        assert khmer.forward_hash(GG, 12) == 11184810

        hi = self.hi
        hi.consume(GG)

        try:
            hi.get(float(GGhash))
            assert "the previous statement should fail"
        except ValueError as err:
            print(str(err))

    def test_collision_1(self):

        GG = 'G' * 12                   # forward_hash: 11184810
        GGhash = khmer.forward_hash(GG, 12)
        assert khmer.forward_hash(GG, 12) == 11184810

        collision_1 = 'AAACGTATGACT'
        assert khmer.forward_hash(collision_1, 12) == 184777

        collision_2 = 'AAATACCGAGCG'
        assert khmer.forward_hash(collision_2, 12) == 76603

        # note, hash(GG) % 1000003 == hash(collision_1)
        # note, hash(GG) % 1009837 == hash(collision_2)

        hi = self.hi
        hi.consume(GG)
        hi.consume(collision_1)

        assert hi.get(GG) == 1
        assert hi.get(GGhash) == 1

    def test_collision_2(self):

        GG = 'G' * 12                   # forward_hash: 11184810
        assert khmer.forward_hash(GG, 12) == 11184810

        collision_1 = 'AAACGTATGACT'
        assert khmer.forward_hash(collision_1, 12) == 184777

        collision_2 = 'AAATACCGAGCG'
        assert khmer.forward_hash(collision_2, 12) == 76603

        # hash(GG) % 1000003 == hash(collision_1)
        # hash(GG) % 1009837 == hash(collision_2)

        hi = self.hi
        hi.consume(GG)
        hi.consume(collision_2)

        assert hi.get(GG) == 1

    def test_collision_3(self):

        GG = 'G' * 12                   # forward_hash: 11184810
        assert khmer.forward_hash(GG, 12) == 11184810

        collision_1 = 'AAACGTATGACT'
        assert khmer.forward_hash(collision_1, 12) == 184777

        collision_2 = 'AAATACCGAGCG'
        assert khmer.forward_hash(collision_2, 12) == 76603

        # hash(GG) % 1000003 == hash(collision_1)
        # hash(GG) % 1009837 == hash(collision_2)

        hi = self.hi
        hi.consume(GG)
        hi.consume(collision_1)
        hi.consume(collision_2)

        assert hi.get(GG) == 2


def test_get_raw_tables():
    ht = khmer.new_counting_hash(20, 1e5, 4)
    tables = ht.get_raw_tables()

    for size, table in zip(ht.hashsizes(), tables):
        assert isinstance(table, memoryview)
        assert size == len(table)


def test_get_raw_tables_view():
    ht = khmer.new_counting_hash(20, 1e5, 4)
    tables = ht.get_raw_tables()
    for tab in tables:
        assert sum(tab.tolist()) == 0
    ht.consume('AAAATTTTCCCCGGGGAAAA')
    for tab in tables:
        assert sum(tab.tolist()) == 1


@attr('huge')
def test_toobig():
    try:
        ct = khmer.new_counting_hash(30, 1e13, 1)
        assert 0, "this should fail"
    except MemoryError as err:
        print(str(err))


def test_3_tables():
    x = list(PRIMES_1m)
    x.append(1000005)

    hi = khmer.CountingHash(12, x)

    GG = 'G' * 12                   # forward_hash: 11184810
    assert khmer.forward_hash(GG, 12) == 11184810

    collision_1 = 'AAACGTATGACT'
    assert khmer.forward_hash(collision_1, 12) == 184777

    collision_2 = 'AAATACCGAGCG'
    assert khmer.forward_hash(collision_2, 12) == 76603

    collision_3 = 'AAACGTATCGAG'
    assert khmer.forward_hash(collision_3, 12) == 184755

    # hash(GG) % 1000003 == hash(collision_1)
    # hash(GG) % 1009837 == hash(collision_2)
    # hash(GG) % 1000005 == hash(collision_3)
    hi.consume(GG)
    assert hi.get(GG) == 1

    hi.consume(collision_1)
    assert hi.get(GG) == 1

    hi.consume(collision_2)
    assert hi.get(GG) == 1

    hi.consume(collision_3)
    assert hi.get(GG) == 2


def test_simple_median():
    hi = khmer.new_counting_hash(6, 1e6, 2)

    hi.consume("AAAAAA")
    (median, average, stddev) = hi.get_median_count("AAAAAA")
    print(median, average, stddev)
    assert median == 1
    assert average == 1.0
    assert stddev == 0.0

    hi.consume("AAAAAA")
    (median, average, stddev) = hi.get_median_count("AAAAAA")
    print(median, average, stddev)
    assert median == 2
    assert average == 2.0
    assert stddev == 0.0

    hi.consume("AAAAAT")
    (median, average, stddev) = hi.get_median_count("AAAAAAT")
    print(median, average, stddev)
    assert median == 2
    assert average == 1.5
    assert int(stddev * 100) == 50        # .5

    hi.consume("AAAAAT")
    (median, average, stddev) = hi.get_median_count("AAAAAAT")
    print(median, average, stddev)
    assert median == 2
    assert average == 2.0
    assert stddev == 0.0

    hi.consume("AAAAAT")
    (median, average, stddev) = hi.get_median_count("AAAAAAT")
    print(median, average, stddev)
    assert median == 3
    assert average == 2.5
    assert int(stddev * 100) == 50        # .5


def test_median_too_short():
    hi = khmer.new_counting_hash(6, 1e6, 2)

    hi.consume("AAAAAA")
    try:
        hi.get_median_count("A")
        assert 0, "this should fail"
    except ValueError:
        pass


def test_median_at_least():
    hi = khmer.new_counting_hash(6, 1e6, 2)

    hi.consume("AAAAAA")
    assert hi.median_at_least("AAAAAA", 1)
    assert hi.median_at_least("AAAAAA", 2) is False

    hi.consume("AAAAAA")
    assert hi.median_at_least("AAAAAA", 2)
    assert hi.median_at_least("AAAAAA", 3) is False

    hi.consume("AAAAAA")
    assert hi.median_at_least("AAAAAA", 3)
    assert hi.median_at_least("AAAAAA", 4) is False

    hi.consume("AAAAAA")
    assert hi.median_at_least("AAAAAA", 4)
    assert hi.median_at_least("AAAAAA", 5) is False

    hi.consume("AAAAAA")
    assert hi.median_at_least("AAAAAA", 5)
    assert hi.median_at_least("AAAAAA", 6) is False


def test_median_at_least_single_gt():
    K = 20
    hi = khmer.new_counting_hash(K, 1e6, 2)

    kmers = ['ATCGATCGATCGATCGATCG',
             'GTACGTACGTACGTACGTAC',
             'TTAGTTAGTTAGTTAGTTAG']

    for kmer in kmers:
        hi.consume(kmer)
        assert hi.median_at_least(kmer, 1) is True


def test_median_at_least_single_lt():
    K = 20
    hi = khmer.new_counting_hash(K, 1e6, 2)

    kmers = ['ATCGATCGATCGATCGATCG',
             'GTACGTACGTACGTACGTAC',
             'TTAGTTAGTTAGTTAGTTAG']

    for kmer in kmers:
        hi.consume(kmer)
        assert hi.median_at_least(kmer, 2) is False


def test_median_at_least_odd_gt():
    # test w/odd number of k-mers
    K = 20
    hi = khmer.new_counting_hash(K, 1e6, 2)

    seqs = ['ATCGATCGATCGATCGATCGCC',
            'GTACGTACGTACGTACGTACCC',
            'TTAGTTAGTTAGTTAGTTAGCC']

    for seq in seqs:
        hi.consume(seq)
        assert hi.median_at_least(seq, 1) is True


def test_median_at_least_odd_lt():
    K = 20
    hi = khmer.new_counting_hash(K, 1e6, 2)

    seqs = ['ATCGATCGATCGATCGATCGCC',
            'GTACGTACGTACGTACGTACCC',
            'TTAGTTAGTTAGTTAGTTAGCC']

    for seq in seqs:
        hi.consume(seq)
        assert hi.median_at_least(seq, 2) is False


# Test median with even number of k-mers
def test_median_at_least_even_gt():
    K = 20
    hi = khmer.new_counting_hash(K, 1e6, 2)

    seqs = ['ATCGATCGATCGATCGATCGCCC',
            'GTACGTACGTACGTACGTACCCC',
            'TTAGTTAGTTAGTTAGTTAGCCC']

    for seq in seqs:
        hi.consume(seq)
        assert hi.median_at_least(seq, 1) is True


def test_median_at_least_even_lt():
    K = 20
    hi = khmer.new_counting_hash(K, 1e6, 2)

    seqs = ['ATCGATCGATCGATCGATCGCCC',
            'GTACGTACGTACGTACGTACCCC',
            'TTAGTTAGTTAGTTAGTTAGCCC']

    for seq in seqs:
        hi.consume(seq)
        assert hi.median_at_least(seq, 2) is False


def test_median_at_least_comp():
    K = 20
    C = 4
    hi = khmer.new_counting_hash(K, 1e6, 2)

    seqs = ['ATCGATCGATCGATCGATCGCCC',
            'GTACGTACGTACGTACGTACCCC',
            'TTAGTTAGTTAGTTAGTTAGCCC']

    for seq in seqs:
        hi.consume(seq)
        hi.consume(seq)
        hi.consume(seq)

        med, _, _ = hi.get_median_count(seq)
        assert hi.median_at_least(seq, C) is (med >= C)


def test_median_at_least_exception():
    ht = khmer.new_counting_hash(20, 1e6, 2)
    try:
        ht.median_at_least('ATGGCTGATCGAT', 1)
        assert 0, "should have thrown ValueError"
    except ValueError as e:
        pass


def test_simple_kadian():
    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    assert hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG") == 1

    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    hi.consume("ACTGCTATCTCTAGAcCTATG")
    #           ---------------^
    x = hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG")
    assert x == 2, x

    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    hi.consume("ACTGCTATCTCTAGAcCTATG")
    #           ---------------^---^
    x = hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG")
    assert x == 2

    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    hi.consume("ACTGCTATCTCTAGtGCTAcG")
    #           --------------^^---^
    x = hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG")
    assert x == 1, x


def test_simple_kadian_2():
    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    assert hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG") == 1

    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    # hi.consume("ACaGCTATCTCTAGAGCTATG")
    hi.consume("ACAGCTATCTCTAGAGCTATG")
    #           --^
    x = hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG")
    assert x == 2, x

    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    # hi.consume("ACaGCTATCTCTAGAcCTATG")
    hi.consume("ACAGCTATCTCTAGACCTATG")
    #           --^          --^
    x = hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG")
    assert x == 1, x

    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    # hi.consume("ACTGCTATCgCTAGAGCTATG")
    hi.consume("ACTGCTATCGCTAGAGCTATG")
    #                  --^
    x = hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG")
    assert x == 2, x


def test_2_kadian():
    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    assert hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG", 2) == 1

    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    # hi.consume("ACTGCTATCTCTAGAcCTATG")
    hi.consume("ACTGCTATCTCTAGACCTATG")
    #           ---------------^
    x = hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG", 2)
    assert x == 2, x

    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    # hi.consume("ACTGCTATCTCTAGAcCTAtG")
    hi.consume("ACTGCTATCTCTAGACCTATG")
    #           ---------------^---^
    assert hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG", 2) == 2

    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    # hi.consume("ACTGCTATCTCTACtcCTAtG")
    hi.consume("ACTGCTATCTCTACTCCTATG")
    #           --------------^^---^
    x = hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG", 2)
    assert x == 2, x

    hi = khmer.new_counting_hash(6, 1e6, 2)
    hi.consume("ACTGCTATCTCTAGAGCTATG")
    # hi.consume("ACTGCTgTCTCTACtcCTAtG")
    hi.consume("ACTGCTGTCTCTACTCCTATG")
    #           ------^-------^^---^
    x = hi.get_kadian_count("ACTGCTATCTCTAGAGCTATG", 2)
    assert x == 1, x


def test_get_kmer_counts_too_short():
    hi = khmer.new_counting_hash(6, 1e6, 2)

    hi.consume("AAAAAA")
    counts = hi.get_kmer_counts("A")
    assert len(counts) == 0


def test_get_kmer_counts_too_short():
    hi = khmer.new_counting_hash(6, 1e6, 2)

    hi.consume("AAAAAA")
    counts = hi.get_kmer_counts("A")
    assert len(counts) == 0


def test_get_kmers_too_short():
    hi = khmer.new_counting_hash(6, 1e6, 2)

    hi.consume("AAAAAA")
    kmers = hi.get_kmers("A")
    assert len(kmers) == 0


def test_get_kmer_counts():
    hi = khmer.new_counting_hash(6, 1e6, 2)

    hi.consume("AAAAAA")
    counts = hi.get_kmer_counts("AAAAAA")
    print(counts)
    assert len(counts) == 1
    assert counts[0] == 1

    hi.consume("AAAAAA")
    counts = hi.get_kmer_counts("AAAAAA")
    print(counts)
    assert len(counts) == 1
    assert counts[0] == 2

    hi.consume("AAAAAT")
    counts = hi.get_kmer_counts("AAAAAAT")
    print(counts)
    assert len(counts) == 2
    assert counts[0] == 2
    assert counts[1] == 1

    hi.consume("AAAAAT")
    counts = hi.get_kmer_counts("AAAAAAT")
    print(counts)
    assert len(counts) == 2
    assert counts[0] == 2
    assert counts[1] == 2

    hi.consume("AAAAAT")
    counts = hi.get_kmer_counts("AAAAAAT")
    print(counts)
    assert len(counts) == 2
    assert counts[0] == 2
    assert counts[1] == 3


def test_get_kmer_hashes():
    hi = khmer.new_counting_hash(6, 1e6, 2)

    hi.consume("AAAAAA")
    hashes = hi.get_kmer_hashes("AAAAAA")
    print(hashes)
    assert len(hashes) == 1
    assert hi.get(hashes[0]) == 1

    hi.consume("AAAAAA")
    hashes = hi.get_kmer_hashes("AAAAAA")
    print(hashes)
    assert len(hashes) == 1
    assert hi.get(hashes[0]) == 2

    hi.consume("AAAAAT")
    hashes = hi.get_kmer_hashes("AAAAAAT")
    print(hashes)
    assert len(hashes) == 2
    assert hi.get(hashes[0]) == 2
    assert hi.get(hashes[1]) == 1

    hi.consume("AAAAAT")
    hashes = hi.get_kmer_hashes("AAAAAAT")
    print(hashes)
    assert len(hashes) == 2
    assert hi.get(hashes[0]) == 2
    assert hi.get(hashes[1]) == 2

    hi.consume("AAAAAT")
    hashes = hi.get_kmer_hashes("AAAAAAT")
    print(hashes)
    assert len(hashes) == 2
    assert hi.get(hashes[0]) == 2
    assert hi.get(hashes[1]) == 3


def test_get_kmers():
    hi = khmer.new_counting_hash(6, 1e6, 2)

    kmers = hi.get_kmers("AAAAAA")
    assert kmers == ["AAAAAA"]

    kmers = hi.get_kmers("AAAAAAT")
    assert kmers == ["AAAAAA", "AAAAAT"]


@attr("huge")
def test_save_load_large():
    def do_test(ctfile):
        inpath = utils.get_test_data('random-20-a.fa')
        savepath = utils.get_temp_filename(ctfile)

        sizes = khmer.get_n_primes_above_x(1, 2**31)

        orig = khmer.CountingHash(12, sizes)
        orig.consume_fasta(inpath)
        orig.save(savepath)

        loaded = khmer.load_counting_hash(savepath)

        orig_count = orig.n_occupied()
        loaded_count = loaded.n_occupied()
        assert orig_count == 3966, orig_count
        assert loaded_count == orig_count, loaded_count

    for ctfile in ['temp.ct.gz', 'temp.ct']:
        do_test(ctfile)


def test_save_load():
    inpath = utils.get_test_data('random-20-a.fa')
    savepath = utils.get_temp_filename('tempcountingsave0.ht')

    sizes = list(PRIMES_1m)
    sizes.append(1000005)

    hi = khmer.CountingHash(12, sizes)
    hi.consume_fasta(inpath)
    hi.save(savepath)

    ht = khmer.CountingHash(12, sizes)
    ht.load(savepath)

    tracking = khmer._Hashbits(12, sizes)
    x = hi.abundance_distribution(inpath, tracking)

    tracking = khmer._Hashbits(12, sizes)
    y = ht.abundance_distribution(inpath, tracking)

    assert sum(x) == 3966, sum(x)
    assert x == y, (x, y)


def test_load_truncated():
    inpath = utils.get_test_data('random-20-a.fa')
    savepath = utils.get_temp_filename('save.ht')
    truncpath = utils.get_temp_filename('trunc.ht')

    sizes = khmer.get_n_primes_near_x(3, 200)

    hi = khmer.CountingHash(12, sizes)
    hi.consume_fasta(inpath)
    hi.save(savepath)

    data = open(savepath, 'rb').read()
    for i in range(len(data)):
        fp = open(truncpath, 'wb')
        fp.write(data[:i])
        fp.close()

        try:
            ht = khmer.load_counting_hash(truncpath)
            assert 0, "this should not be reached!"
        except IOError as err:
            print(str(err))


def test_load_gz():
    inpath = utils.get_test_data('random-20-a.fa')

    savepath = utils.get_temp_filename('tempcountingsave1.ht')
    loadpath = utils.get_temp_filename('tempcountingsave1.ht.gz')

    sizes = list(PRIMES_1m)
    sizes.append(1000005)

    # save uncompressed hashtable.
    hi = khmer.CountingHash(12, sizes)
    hi.consume_fasta(inpath)
    hi.save(savepath)

    # compress.
    in_file = open(savepath, 'rb')
    out_file = gzip.open(loadpath, 'wb')
    out_file.writelines(in_file)
    out_file.close()
    in_file.close()

    # load compressed hashtable.
    ht = khmer.CountingHash(12, sizes)
    ht.load(loadpath)

    tracking = khmer._Hashbits(12, sizes)
    x = hi.abundance_distribution(inpath, tracking)

    tracking = khmer._Hashbits(12, sizes)
    y = ht.abundance_distribution(inpath, tracking)

    assert sum(x) == 3966, sum(x)
    assert x == y, (x, y)


def test_save_load_gz():
    inpath = utils.get_test_data('random-20-a.fa')
    savepath = utils.get_temp_filename('tempcountingsave2.ht.gz')

    sizes = list(PRIMES_1m)
    sizes.append(1000005)

    hi = khmer.CountingHash(12, sizes)
    hi.consume_fasta(inpath)
    hi.save(savepath)

    ht = khmer.CountingHash(12, sizes)
    ht.load(savepath)

    tracking = khmer._Hashbits(12, sizes)
    x = hi.abundance_distribution(inpath, tracking)

    tracking = khmer._Hashbits(12, sizes)
    y = ht.abundance_distribution(inpath, tracking)

    assert sum(x) == 3966, sum(x)
    assert x == y, (x, y)


def test_trim_full():
    hi = khmer.new_counting_hash(6, 1e6, 2)

    hi.consume(DNA)
    hi.consume(DNA)

    seq, pos = hi.trim_on_abundance(DNA, 2)
    assert DNA == seq, seq


def test_trim_short():
    hi = khmer.new_counting_hash(6, 1e6, 2)

    hi.consume(DNA)
    hi.consume(DNA[:50])

    seq, pos = hi.trim_on_abundance(DNA, 2)
    assert DNA[:50] == seq, (seq, pos)
    assert hi.get(seq[-6:]) == 2
    assert hi.get(DNA[:51][-6:]) == 1


def test_find_spectral_error_positions_1():
    hi = khmer.new_counting_hash(8, 1e6, 2)

    hi.consume(DNA)
    hi.consume(DNA[:30])

    for n in range(len(DNA) - 8 + 1):
        print(n, hi.get(DNA[n:n + 8]))

    posns = hi.find_spectral_error_positions(DNA, 1)
    assert posns == [30], posns


def test_find_spectral_error_positions_2():
    hi = khmer.new_counting_hash(8, 1e6, 2)

    hi.consume(DNA)
    hi.consume(DNA)

    posns = hi.find_spectral_error_positions(DNA, 2)
    assert posns == [], posns


def test_find_spectral_error_positions_6():
    hi = khmer.new_counting_hash(8, 1e6, 2)

    hi.consume(DNA)
    hi.consume(DNA[1:])

    for n in range(len(DNA) - 8 + 1):
        print(n, hi.get(DNA[n:n + 8]))

    posns = hi.find_spectral_error_positions(DNA, 1)
    assert posns == [0], posns


def test_find_spectral_error_positions_4():
    hi = khmer.new_counting_hash(8, 1e6, 2)

    hi.consume(DNA)

    posns = hi.find_spectral_error_positions(DNA, 2)
    assert posns == [], posns


def test_find_spectral_error_positions_5():
    hi = khmer.new_counting_hash(8, 1e6, 2)

    hi.consume(DNA)
    hi.consume(DNA[:10])
    hi.consume(DNA[11:])

    posns = hi.find_spectral_error_positions(DNA, 1)
    assert posns == [10], posns


def test_find_spectral_error_positions_6():
    K = 8
    hi = khmer.new_counting_hash(K, 1e6, 2)

    hi.consume(DNA)
    hi.consume(DNA[K:])

    for n in range(len(DNA) - 8 + 1):
        print(n, hi.get(DNA[n:n + 8]))

    posns = hi.find_spectral_error_positions(DNA, 1)
    assert posns == [7], posns


def test_find_spectral_error_positions_err():
    hi = khmer.new_counting_hash(8, 1e6, 2)

    try:
        posns = hi.find_spectral_error_positions(DNA[:6], 1)
        assert 0, "should raise ValueError; too short"
    except ValueError:
        pass

    try:
        posns = hi.find_spectral_error_positions("ACGTACGN", 1)
        assert 0, "should raise ValueError; contains N"
    except ValueError:
        pass


def test_maxcount():
    # hashtable should saturate at some point so as not to overflow counter
    kh = khmer.new_counting_hash(4, 4 ** 4, 4)
    kh.set_use_bigcount(False)

    last_count = None
    for i in range(0, 1000):
        kh.count('AAAA')
        c = kh.get('AAAA')

        if c == last_count:
            break
        last_count = c

    assert c != 1000, "should not be able to count to 1000: %d" % c
    assert c == MAX_COUNT, c       # this will depend on HashcountType...


def test_maxcount_with_bigcount():
    # hashtable should not saturate, if use_bigcount is set.
    kh = khmer.new_counting_hash(4, 4 ** 4, 4)
    kh.set_use_bigcount(True)

    last_count = None
    for i in range(0, 1000):
        kh.count('AAAA')
        c = kh.get('AAAA')

        if c == last_count:
            break
        last_count = c

    assert c == 1000, "should be able to count to 1000: %d" % c
    assert c != MAX_COUNT, c


def test_maxcount_with_bigcount_save():
    # hashtable should not saturate, if use_bigcount is set.
    kh = khmer.new_counting_hash(4, 4 ** 4, 4)
    kh.set_use_bigcount(True)

    for i in range(0, 1000):
        kh.count('AAAA')
        c = kh.get('AAAA')

    savepath = utils.get_temp_filename('tempcountingsave.ht')
    kh.save(savepath)

    kh = khmer.new_counting_hash(1, 1, 1)
    kh.load(savepath)

    c = kh.get('AAAA')
    assert c == 1000, "should be able to count to 1000: %d" % c
    assert c != MAX_COUNT, c


def test_bigcount_save():
    # hashtable should not saturate, if use_bigcount is set.
    kh = khmer.new_counting_hash(4, 4 ** 4, 4)
    kh.set_use_bigcount(True)

    savepath = utils.get_temp_filename('tempcountingsave.ht')
    kh.save(savepath)

    kh = khmer.new_counting_hash(1, 1, 1)
    kh.load(savepath)

    # set_use_bigcount should still be True after load (i.e. should be saved)

    assert kh.get('AAAA') == 0

    for i in range(0, 1000):
        kh.count('AAAA')
        kh.get('AAAA')

    assert kh.get('AAAA') == 1000


def test_nobigcount_save():
    kh = khmer.new_counting_hash(4, 4 ** 4, 4)
    # kh.set_use_bigcount(False) <-- this is the default

    savepath = utils.get_temp_filename('tempcountingsave.ht')
    kh.save(savepath)

    kh = khmer.new_counting_hash(1, 1, 1)
    kh.load(savepath)

    # set_use_bigcount should still be False after load (i.e. should be saved)

    assert kh.get('AAAA') == 0

    for i in range(0, 1000):
        kh.count('AAAA')
        kh.get('AAAA')

    assert kh.get('AAAA') == MAX_COUNT


def test_bigcount_abund_dist():
    kh = khmer.new_counting_hash(18, 1e2, 4)
    tracking = khmer.new_hashbits(18, 1e2, 4)
    kh.set_use_bigcount(True)

    seqpath = utils.get_test_data('test-abund-read-2.fa')

    kh.consume_fasta(seqpath)

    dist = kh.abundance_distribution(seqpath, tracking)
    print(kh.get('GGTTGACGGGGCTCAGGG'))

    pdist = [(i, dist[i]) for i in range(len(dist)) if dist[i]]
    assert dist[1001] == 1, pdist


def test_bigcount_abund_dist_2():
    kh = khmer.new_counting_hash(18, 1e7, 4)
    tracking = khmer.new_hashbits(18, 1e7, 4)
    kh.set_use_bigcount(True)

    seqpath = utils.get_test_data('test-abund-read.fa')

    kh.consume_fasta(seqpath)
    for i in range(1000):
        kh.count('GGTTGACGGGGCTCAGGG')

    dist = kh.abundance_distribution(seqpath, tracking)
    print(kh.get('GGTTGACGGGGCTCAGGG'))

    pdist = [(i, dist[i]) for i in range(len(dist)) if dist[i]]
    assert dist[1001] == 1, pdist


def test_bigcount_overflow():
    kh = khmer.new_counting_hash(18, 1e7, 4)
    kh.set_use_bigcount(True)

    for i in range(0, 70000):
        kh.count('GGTTGACGGGGCTCAGGG')

    assert kh.get('GGTTGACGGGGCTCAGGG') == MAX_BIGCOUNT


def test_get_ksize():
    kh = khmer.new_counting_hash(22, 1, 1)
    assert kh.ksize() == 22


def test_get_hashsizes():
    kh = khmer.new_counting_hash(22, 100, 4)
    assert kh.hashsizes() == [101, 103, 107, 109], kh.hashsizes()

# def test_collect_high_abundance_kmers():
#    seqpath = utils.get_test_data('test-abund-read-2.fa')
#
#    kh = khmer.new_counting_hash(18, 1e6, 4)
#    hb = kh.collect_high_abundance_kmers(seqpath, 2, 4)


#


def test_load_notexist_should_fail():
    savepath = utils.get_temp_filename('tempcountingsave0.ht')

    hi = khmer.new_counting_hash(12, 1000)
    try:
        hi.load(savepath)
        assert 0, "load should fail"
    except IOError as e:
        print(str(e))


def test_load_truncated_should_fail():
    inpath = utils.get_test_data('random-20-a.fa')
    savepath = utils.get_temp_filename('tempcountingsave0.ht')

    hi = khmer.new_counting_hash(12, 1000)
    hi.consume_fasta(inpath)
    hi.save(savepath)

    fp = open(savepath, 'rb')
    data = fp.read()
    fp.close()

    fp = open(savepath, 'wb')
    fp.write(data[:1000])
    fp.close()

    hi = khmer.new_counting_hash(12, 1)
    try:
        hi.load(savepath)
        assert 0, "load should fail"
    except IOError as e:
        print(str(e))


def test_load_gz_notexist_should_fail():
    savepath = utils.get_temp_filename('tempcountingsave0.ht.gz')

    hi = khmer.new_counting_hash(12, 1000)
    try:
        hi.load(savepath)
        assert 0, "load should fail"
    except IOError as e:
        print(str(e))


def test_load_gz_truncated_should_fail():
    inpath = utils.get_test_data('random-20-a.fa')
    savepath = utils.get_temp_filename('tempcountingsave0.ht.gz')

    hi = khmer.new_counting_hash(12, 1000)
    hi.consume_fasta(inpath)
    hi.save(savepath)

    fp = open(savepath, 'rb')
    data = fp.read()
    fp.close()

    fp = open(savepath, 'wb')
    fp.write(data[:1000])
    fp.close()

    hi = khmer.new_counting_hash(12, 1)
    try:
        hi.load(savepath)
        assert 0, "load should fail"
    except IOError as e:
        print(str(e))


def test_counting_file_version_check():
    ht = khmer.new_counting_hash(12, 1, 1)

    inpath = utils.get_test_data('badversion-k12.ct')

    try:
        ht.load(inpath)
        assert 0, "this should fail"
    except IOError as e:
        print(str(e))


def test_counting_gz_file_version_check():
    ht = khmer.new_counting_hash(12, 1, 1)

    inpath = utils.get_test_data('badversion-k12.ct.gz')

    try:
        ht.load(inpath)
        assert 0, "this should fail"
    except IOError as e:
        print(str(e))


def test_counting_file_type_check():
    inpath = utils.get_test_data('goodversion-k12.ht')

    kh = khmer.new_counting_hash(12, 1, 1)

    try:
        kh.load(inpath)
        assert 0, "this should fail"
    except IOError as e:
        print(str(e))


def test_counting_gz_file_type_check():
    ht = khmer.new_hashbits(12, 1, 1)

    inpath = utils.get_test_data('goodversion-k12.ht.gz')

    kh = khmer.new_counting_hash(12, 1, 1)

    try:
        kh.load(inpath)
        assert 0, "this should fail"
    except IOError as e:
        print(str(e))


def test_counting_bad_primes_list():
    try:
        ht = khmer.CountingHash(12, ["a", "b", "c"], 1)
        assert 0, "bad list of primes should fail"
    except TypeError as e:
        print(str(e))


def test_bad_use_bigcount():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    countingtable.set_use_bigcount(True)
    assert countingtable.get_use_bigcount()
    try:
        countingtable.get_use_bigcount(True)
        assert 0, "this should fail"
    except TypeError as err:
        print(str(err))


def test_consume_absentfasta():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.consume_fasta("absent_file.fa")
        assert 0, "This should fail"
    except IOError as err:
        print(str(err))


def test_consume_absentfasta_with_reads_parser():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.consume_fasta_with_reads_parser()
        assert 0, "this should fail"
    except TypeError as err:
        print(str(err))
    try:
        readparser = ReadParser(utils.get_test_data('empty-file'))
        countingtable.consume_fasta_with_reads_parser(readparser)
        assert 0, "this should fail"
    except IOError as err:
        print(str(err))
    except ValueError as err:
        print(str(err))


def test_badconsume():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.consume()
        assert 0, "this should fail"
    except TypeError as err:
        print(str(err))
    try:
        countingtable.consume("AAA")
        assert 0, "this should fail"
    except ValueError as err:
        print(str(err))


def test_get_badmin_count():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.get_min_count()
        assert 0, "this should fail"
    except TypeError as err:
        print(str(err))
    try:
        countingtable.get_min_count("AAA")
        assert 0, "this should fail"
    except ValueError as err:
        print(str(err))


def test_get_badmax_count():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.get_max_count()
        assert 0, "this should fail"
    except TypeError as err:
        print(str(err))
    try:
        countingtable.get_max_count("AAA")
        assert 0, "this should fail"
    except ValueError as err:
        print(str(err))


def test_get_badmedian_count():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.get_median_count()
        assert 0, "this should fail"
    except TypeError as err:
        print(str(err))
    try:
        countingtable.get_median_count("AAA")
        assert 0, "this should fail"
    except ValueError as err:
        print(str(err))


def test_get_badkadian_count():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.get_kadian_count()
        assert 0, "this should fail"
    except TypeError as err:
        print(str(err))
    try:
        countingtable.get_kadian_count("AAA")
        assert 0, "this should fail"
    except ValueError as err:
        print(str(err))


def test_badget():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.get()
        assert 0, "this should fail"
    except TypeError as err:
        print(str(err))


def test_badget_2():
    countingtable = khmer.new_counting_hash(6, 1e6)

    countingtable.consume(DNA)

    assert countingtable.get("AGCTTT") == 1

    assert countingtable.get("GATGAG") == 0

    try:
        countingtable.get("AGCTT")
        assert 0, "this should fail"
    except ValueError as err:
        print(str(err))


def test_badtrim():
    countingtable = khmer.new_counting_hash(6, 1e6, 2)

    countingtable.consume(DNA)
    try:
        countingtable.trim_on_abundance()
        assert 0, "this should fail"
    except TypeError as err:
        print(str(err))
    countingtable.trim_on_abundance("AAAAAA", 1)


def test_badfasta_count_kmers_by_position():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.fasta_count_kmers_by_position()
    except TypeError as err:
        print(str(err))

    filename = utils.get_test_data("test-short.fa")
    try:
        countingtable.fasta_count_kmers_by_position(filename, -1, 0)
        assert 0, "this should fail"
    except ValueError as err:
        print(str(err))
    try:
        countingtable.fasta_count_kmers_by_position(filename, 0, -1)
        assert 0, "this should fail"
    except ValueError as err:
        print(str(err))


def test_badload():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.load()
        assert 0, "this should fail"
    except TypeError as err:
        print(str(err))


def test_badsave():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.save()
        assert 0, "this should fail"
    except TypeError as err:
        print(str(err))


def test_badksize():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.ksize(True)
        assert 0, "this should fail"
    except TypeError as err:
        print(str(err))


def test_badhashsizes():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.hashsizes(True)
        assert 0, "this should fail"
    except TypeError as err:
        print(str(err))


def test_badconsume_and_tag():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.consume_and_tag()
        assert 0, "this should fail"
    except TypeError as err:
        print(str(err))


def test_consume_fasta_and_tag():
    countingtable = khmer.new_counting_hash(4, 4 ** 4, 4)
    try:
        countingtable.consume_fasta_and_tag()
        assert 0, "this should fail"
    except TypeError as err:
        print(str(err))
    countingtable.consume_fasta_and_tag(utils.get_test_data("test-graph2.fa"))


def test_consume_and_retrieve_tags_1():
    ct = khmer.new_counting_hash(4, 4 ** 4, 4)

    # first, for each sequence, build tags.
    for record in screed.open(utils.get_test_data('test-graph2.fa')):
        ct.consume_and_tag(record.sequence)

    # check that all the tags in sequences are retrieved by iterating
    # across the sequence and retrieving by neighborhood.

    ss = set()
    tt = set()
    for record in screed.open(utils.get_test_data('test-graph2.fa')):
        for p, tag in ct.get_tags_and_positions(record.sequence):
            ss.add(tag)

        for start in range(len(record.sequence) - 3):
            kmer = record.sequence[start:start + 4]
            tt.update(ct.find_all_tags_list(kmer))

    assert ss == tt


def test_consume_and_retrieve_tags_empty():
    ct = khmer.new_counting_hash(4, 4 ** 4, 4)

    # load each sequence but do not build tags - everything should be empty.
    for record in screed.open(utils.get_test_data('test-graph2.fa')):
        ct.consume(record.sequence)

    # check that all the tags in sequences are retrieved by iterating
    # across the sequence and retrieving by neighborhood.

    ss = set()
    tt = set()
    for record in screed.open(utils.get_test_data('test-graph2.fa')):
        for p, tag in ct.get_tags_and_positions(record.sequence):
            ss.add(tag)

        for start in range(len(record.sequence) - 3):
            kmer = record.sequence[start:start + 4]
            tt.update(ct.find_all_tags_list(kmer))

    assert not ss
    assert not tt


def test_find_all_tags_list_error():
    ct = khmer.new_counting_hash(4, 4 ** 4, 4)

    # load each sequence but do not build tags - everything should be empty.
    for record in screed.open(utils.get_test_data('test-graph2.fa')):
        ct.consume(record.sequence)

    try:
        ct.find_all_tags_list("ATA")
        assert False, "a ValueError should be raised for incorrect k-mer size"
    except ValueError:
        pass

    try:
        ct.find_all_tags_list("ATAGA")
        assert False, "a ValueError should be raised for incorrect k-mer size"
    except ValueError:
        pass


def test_abund_dist_gz_bigcount():
    infile = utils.get_temp_filename('test.fa')
    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)
    outfile = utils.get_temp_filename('test_ct.gz')
    script = 'load-into-counting.py'
    htfile = utils.get_temp_filename('test_ct')
    args = ['-x', str(1e7), '-N', str(2), '-k', str(2), htfile, infile]
    utils.runscript(script, args)  # create a bigcount table
    assert os.path.exists(htfile)
    data = open(htfile, 'rb').read()
    f_out = gzip.open(outfile, 'wb')  # compress the created bigcount table
    f_out.write(data)
    f_out.close()
    # load the compressed bigcount table
    counting_hash = khmer.load_counting_hash(outfile)
    hashsizes = counting_hash.hashsizes()
    kmer_size = counting_hash.ksize()
    tracking = khmer._Hashbits(kmer_size, hashsizes)
    abundances = counting_hash.abundance_distribution(infile, tracking)
    # calculate abundance distribution for compressed bigcount table
    flag = False
    # check if abundance is > 255
    # if ok  gzipped bigcount was loaded correctly
    for _, i in enumerate(abundances):
        print(_, i)
        if _ > 255 and i > 0:
            flag = True
            break
    assert flag


def test_counting_load_bigcount():
    count_table = khmer.new_counting_hash(10, 1e5, 4)
    count_table.set_use_bigcount(True)
    for i in range(500):
        print(i, count_table.count('ATATATATAT'))
    count = count_table.get('ATATATATAT')
    assert count == 500
