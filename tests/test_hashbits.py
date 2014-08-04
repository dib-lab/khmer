#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=missing-docstring,protected-access
import khmer
from khmer import ReadParser

from screed.fasta import fasta_iter
import screed

import khmer_tst_utils as utils
from nose.plugins.attrib import attr


def teardown():
    utils.cleanup()


def test__get_set_tag_density():
    ht = khmer.new_hashbits(32, 1, 1)

    orig = ht._get_tag_density()
    assert orig != 2
    ht._set_tag_density(2)
    assert ht._get_tag_density() == 2


def test_n_occupied_1():
    filename = utils.get_test_data('random-20-a.fa')

    K = 20  # size of kmer
    HT_SIZE = 100000  # size of hashtable
    N_HT = 1  # number of hashtables

    # test modified c++ n_occupied code
    ht1 = khmer.new_hashbits(K, HT_SIZE, N_HT)

    for n, record in enumerate(fasta_iter(open(filename))):
        ht1.consume(record['sequence'])

    # this number calculated independently
    assert ht1.n_occupied() == 3877


def test_bloom_python_1():
    # test python code to count unique kmers using bloom filter
    filename = utils.get_test_data('random-20-a.fa')

    K = 20  # size of kmer
    HT_SIZE = 100000  # size of hashtable
    N_HT = 3  # number of hashtables

    ht2 = khmer.new_hashbits(K, HT_SIZE, N_HT)

    n_unique = 0
    for n, record in enumerate(fasta_iter(open(filename))):
        sequence = record['sequence']
        seq_len = len(sequence)
        for n in range(0, seq_len + 1 - K):
            kmer = sequence[n:n + K]
            if (not ht2.get(kmer)):
                n_unique += 1
            ht2.count(kmer)

    assert n_unique == 3960
    assert ht2.n_occupied() == 3882
    assert ht2.n_unique_kmers() == 3960  # this number equals to n_unique


def test_bloom_c_1():
    # test c++ code to count unique kmers using bloom filter

    filename = utils.get_test_data('random-20-a.fa')

    K = 20  # size of kmer
    HT_SIZE = 100000  # size of hashtable
    N_HT = 3  # number of hashtables

    ht3 = khmer.new_hashbits(K, HT_SIZE, N_HT)

    for n, record in enumerate(fasta_iter(open(filename))):
        ht3.consume(record['sequence'])

    assert ht3.n_occupied() == 3882
    assert ht3.n_unique_kmers() == 3960


def test_n_occupied_2():  # simple one
    K = 4
    HT_SIZE = 10  # use 11
    N_HT = 1

    ht1 = khmer.new_hashbits(K, HT_SIZE, N_HT)
    ht1.count('AAAA')  # 00 00 00 00 = 0
    assert ht1.n_occupied() == 1

    ht1.count('ACTG')  # 00 10 01 11 =
    assert ht1.n_occupied() == 2

    ht1.count('AACG')  # 00 00 10 11 = 11  # collision 1

    assert ht1.n_occupied() == 2
    ht1.count('AGAC')   # 00  11 00 10 # collision 2
    assert ht1.n_occupied() == 2


def test_bloom_c_2():  # simple one
    K = 4
    HT_SIZE = 10  # use 11
    N_HT1 = 1    # hashtable size = 11
    N_HT2 = 2    # hashtable size = 11,13

    # use only 1 hashtable, no bloom filter
    ht1 = khmer.new_hashbits(K, HT_SIZE, N_HT1)
    ht1.count('AAAA')  # 00 00 00 00 = 0
    ht1.count('ACTG')  # 00 10 01 11 =
    assert ht1.n_unique_kmers() == 2
    ht1.count('AACG')  # 00 00 10 11 = 11  # collision  with 1st kmer
    assert ht1.n_unique_kmers() == 2
    ht1.count('AGAC')   # 00  11 00 10 # collision  with 2nd kmer
    assert ht1.n_unique_kmers() == 2

    # use two hashtables with 11,13
    ht2 = khmer.new_hashbits(K, HT_SIZE, N_HT2)
    ht2.count('AAAA')  # 00 00 00 00 = 0

    ht2.count('ACTG')  # 00 10 01 11 = 2*16 +4 +3 = 39
    assert ht2.n_unique_kmers() == 2
    ht2.count('AACG')  # 00 00 10 11 = 11  # collision with only 1st kmer
    assert ht2.n_unique_kmers() == 3
    ht2.count('AGAC')   # 00  11 00 10  3*16 +2 = 50
    # collision with both 2nd and 3rd kmers

    assert ht2.n_unique_kmers() == 3


@attr('highmem')
def test_filter_if_present():
    ht = khmer.new_hashbits(32, 2, 2)

    maskfile = utils.get_test_data('filter-test-A.fa')
    inputfile = utils.get_test_data('filter-test-B.fa')
    outfile = utils.get_temp_filename('filter')

    ht.consume_fasta(maskfile)
    ht.filter_if_present(inputfile, outfile)

    records = list(fasta_iter(open(outfile)))
    assert len(records) == 1
    assert records[0]['name'] == '3'


@attr('highmem')
def test_combine_pe():
    inpfile = utils.get_test_data('combine_parts_1.fa')
    ht = khmer.new_hashbits(32, 1, 1)

    ht.consume_partitioned_fasta(inpfile)
    assert ht.count_partitions() == (2, 0)

    s1 = "CATGCAGAAGTTCCGCAACCATACCGTTCAGT"
    pid1 = ht.get_partition_id(s1)

    s2 = "CAAATGTACATGCACTTAAAATCATCCAGCCG"
    pid2 = ht.get_partition_id(s2)

    assert pid1 == 2
    assert pid2 == 80293

    ht.join_partitions(pid1, pid2)

    pid1 = ht.get_partition_id(s1)
    pid2 = ht.get_partition_id(s2)

    assert pid1 == pid2
    assert ht.count_partitions() == (1, 0)


@attr('highmem')
def test_load_partitioned():
    inpfile = utils.get_test_data('combine_parts_1.fa')
    ht = khmer.new_hashbits(32, 1, 1)

    ht.consume_partitioned_fasta(inpfile)
    assert ht.count_partitions() == (2, 0)

    s1 = "CATGCAGAAGTTCCGCAACCATACCGTTCAGT"
    assert ht.get(s1)

    s2 = "CAAATGTACATGCACTTAAAATCATCCAGCCG"
    assert ht.get(s2)

    s3 = "CATGCAGAAGTTCCGCAACCATACCGTTCAGTTCCTGGTGGCTA"[-32:]
    assert ht.get(s3)


@attr('highmem')
def test_count_within_radius_simple():
    inpfile = utils.get_test_data('all-A.fa')
    ht = khmer.new_hashbits(4, 2, 2)

    print ht.consume_fasta(inpfile)
    n = ht.count_kmers_within_radius('AAAA', 1)
    assert n == 1

    n = ht.count_kmers_within_radius('AAAA', 10)
    assert n == 1


@attr('highmem')
def test_count_within_radius_big():
    inpfile = utils.get_test_data('random-20-a.fa')
    ht = khmer.new_hashbits(20, 1e5, 4)

    ht.consume_fasta(inpfile)
    n = ht.count_kmers_within_radius('CGCAGGCTGGATTCTAGAGG', int(1e6))
    assert n == 3960

    ht = khmer.new_hashbits(21, 1e5, 4)
    ht.consume_fasta(inpfile)
    n = ht.count_kmers_within_radius('CGCAGGCTGGATTCTAGAGGC', int(1e6))
    assert n == 39


@attr('highmem')
def test_count_kmer_degree():
    inpfile = utils.get_test_data('all-A.fa')
    ht = khmer.new_hashbits(4, 2, 2)
    ht.consume_fasta(inpfile)

    assert ht.kmer_degree('AAAA') == 2
    assert ht.kmer_degree('AAAT') == 1
    assert ht.kmer_degree('AATA') == 0
    assert ht.kmer_degree('TAAA') == 1


@attr('highmem')
def test_find_radius_for_volume():
    inpfile = utils.get_test_data('all-A.fa')
    ht = khmer.new_hashbits(4, 2, 2)
    ht.consume_fasta(inpfile)

    assert ht.find_radius_for_volume('AAAA', 0, 100) == 0
    assert ht.find_radius_for_volume('AAAA', 1, 100) == 0
    assert ht.find_radius_for_volume('AAAA', 2, 100) == 100


def test_circumference():
    ht = khmer.new_hashbits(4, 10, 2)

    ht.count('ATGC')
    ht.count('GATG')
    ht.count('ATGG')

    x = ht.count_kmers_on_radius('GATG', 1, 200)
    assert x == 2

    ht.count('ATGA')
    x = ht.count_kmers_on_radius('GATG', 1, 200)
    assert x == 3, x

    ht.count('TGAT')
    x = ht.count_kmers_on_radius('GATG', 1, 200)
    assert x == 4, x


def test_save_load_tagset():
    ht = khmer.new_hashbits(32, 1, 1)

    outfile = utils.get_temp_filename('tagset')

    ht.add_tag('A' * 32)
    ht.save_tagset(outfile)

    ht.add_tag('G' * 32)

    ht.load_tagset(outfile)              # implicitly => clear_tags=True
    ht.save_tagset(outfile)

    # if tags have been cleared, then the new tagfile will be larger (34 bytes)
    # else smaller (26 bytes).

    fp = open(outfile, 'rb')
    data = fp.read()
    fp.close()
    assert len(data) == 26, len(data)


def test_save_load_tagset_noclear():
    ht = khmer.new_hashbits(32, 1, 1)

    outfile = utils.get_temp_filename('tagset')

    ht.add_tag('A' * 32)
    ht.save_tagset(outfile)

    ht.add_tag('G' * 32)

    ht.load_tagset(outfile, False)       # set clear_tags => False; zero tags
    ht.save_tagset(outfile)

    # if tags have been cleared, then the new tagfile will be large (34 bytes);
    # else small (26 bytes).

    fp = open(outfile, 'rb')
    data = fp.read()
    fp.close()
    assert len(data) == 34, len(data)


@attr('highmem')
def test_stop_traverse():
    filename = utils.get_test_data('random-20-a.fa')

    K = 20  # size of kmer
    HT_SIZE = 1e4  # size of hashtable
    N_HT = 3  # number of hashtables

    ht = khmer.new_hashbits(K, HT_SIZE, N_HT)

    # without tagging/joining across consume, this breaks into two partition;
    # with, it is one partition.
    ht.add_stop_tag('TTGCATACGTTGAGCCAGCG')

    ht.consume_fasta_and_tag(filename)   # DO NOT join reads across stoptags
    subset = ht.do_subset_partition(0, 0, True)
    ht.merge_subset(subset)

    n, _ = ht.count_partitions()
    assert n == 2, n


@attr('highmem')
def test_tag_across_stoptraverse():
    filename = utils.get_test_data('random-20-a.fa')

    K = 20  # size of kmer
    HT_SIZE = 1e4  # size of hashtable
    N_HT = 3  # number of hashtables

    ht = khmer.new_hashbits(K, HT_SIZE, N_HT)

    # without tagging/joining across consume, this breaks into two partition;
    # with, it is one partition.
    ht.add_stop_tag('CCGAATATATAACAGCGACG')

    ht.consume_fasta_and_tag_with_stoptags(filename)  # DO join reads across

    subset = ht.do_subset_partition(0, 0)
    n, _ = ht.count_partitions()
    assert n == 99                       # reads only connected by traversal...

    n, _ = ht.subset_count_partitions(subset)
    assert n == 2                        # but need main to cross stoptags.

    ht.merge_subset(subset)

    n, _ = ht.count_partitions()         # ta-da!
    assert n == 1, n


@attr('highmem')
def test_notag_across_stoptraverse():
    filename = utils.get_test_data('random-20-a.fa')

    K = 20  # size of kmer
    HT_SIZE = 1e4  # size of hashtable
    N_HT = 3  # number of hashtables

    ht = khmer.new_hashbits(K, HT_SIZE, N_HT)

    # connecting k-mer at the beginning/end of a read: breaks up into two.
    ht.add_stop_tag('TTGCATACGTTGAGCCAGCG')

    ht.consume_fasta_and_tag_with_stoptags(filename)

    subset = ht.do_subset_partition(0, 0)
    ht.merge_subset(subset)

    n, _ = ht.count_partitions()
    assert n == 2, n


def test_find_stoptags():
    ht = khmer.new_hashbits(5, 1, 1)
    ht.add_stop_tag("AAAAA")

    assert ht.identify_stoptags_by_position("AAAAA") == [0]
    assert ht.identify_stoptags_by_position("AAAAAA") == [0, 1]
    assert ht.identify_stoptags_by_position("TTTTT") == [0]
    assert ht.identify_stoptags_by_position("TTTTTT") == [0, 1]


def test_find_stoptags2():
    ht = khmer.new_hashbits(4, 1, 1)
    ht.add_stop_tag("ATGC")

    x = ht.identify_stoptags_by_position("ATGCATGCGCAT")
    assert x == [0, 2, 4, 8], x


def test_get_ksize():
    kh = khmer.new_hashbits(22, 1, 1)
    assert kh.ksize() == 22


def test_get_hashsizes():
    kh = khmer.new_hashbits(22, 100, 4)
    assert kh.hashsizes() == [101, 103, 107, 109], kh.hashsizes()


def test_extract_unique_paths_0():
    kh = khmer.new_hashbits(10, 4, 4)

    x = kh.extract_unique_paths('ATGGAGAGACACAGATAGACAGGAGTGGCGATG', 10, 1)
    assert x == ['ATGGAGAGACACAGATAGACAGGAGTGGCGATG']

    kh.consume('ATGGAGAGACACAGATAGACAGGAGTGGCGATG')
    x = kh.extract_unique_paths('ATGGAGAGACACAGATAGACAGGAGTGGCGATG', 10, 1)
    assert not x


def test_extract_unique_paths_1():
    kh = khmer.new_hashbits(10, 4, 4)

    kh.consume('AGTGGCGATG')
    x = kh.extract_unique_paths('ATGGAGAGACACAGATAGACAGGAGTGGCGATG', 10, 1)
    print x
    assert x == ['ATGGAGAGACACAGATAGACAGGAGTGGCGAT']  # all but the last k-mer


def test_extract_unique_paths_2():
    kh = khmer.new_hashbits(10, 4, 4)

    kh.consume('ATGGAGAGAC')
    x = kh.extract_unique_paths('ATGGAGAGACACAGATAGACAGGAGTGGCGATG', 10, 1)
    print x
    assert x == ['TGGAGAGACACAGATAGACAGGAGTGGCGATG']  # all but the 1st k-mer


def test_extract_unique_paths_3():
    kh = khmer.new_hashbits(10, 4, 4)

    kh.consume('ATGGAGAGAC')
    kh.consume('AGTGGCGATG')
    x = kh.extract_unique_paths('ATGGAGAGACACAGATAGACAGGAGTGGCGATG', 10, 1)
    print x
    # all but the 1st/last k-mer
    assert x == ['TGGAGAGACACAGATAGACAGGAGTGGCGAT']


def test_extract_unique_paths_4():
    kh = khmer.new_hashbits(10, 4, 4)

    kh.consume('ATGGAGAGAC')
    kh.consume('AGTGGCGATG')

    kh.consume('ATAGACAGGA')

    x = kh.extract_unique_paths('ATGGAGAGACACAGATAGACAGGAGTGGCGATG', 10, 1)
    print x
    assert x == ['TGGAGAGACACAGATAGACAGG', 'TAGACAGGAGTGGCGAT']


@attr('highmem')
def test_find_unpart():
    filename = utils.get_test_data('random-20-a.odd.fa')
    filename2 = utils.get_test_data('random-20-a.even.fa')

    K = 20  # size of kmer
    HT_SIZE = 1e4  # size of hashtable
    N_HT = 3  # number of hashtables

    ht = khmer.new_hashbits(K, HT_SIZE, N_HT)
    ht.consume_fasta_and_tag(filename)

    subset = ht.do_subset_partition(0, 0)
    ht.merge_subset(subset)

    n, _ = ht.count_partitions()
    assert n == 49

    ht.find_unpart(filename2, True, False)
    n, _ = ht.count_partitions()
    assert n == 1, n                     # all sequences connect


@attr('highmem')
def test_find_unpart_notraverse():
    filename = utils.get_test_data('random-20-a.odd.fa')
    filename2 = utils.get_test_data('random-20-a.even.fa')

    K = 20  # size of kmer
    HT_SIZE = 1e4  # size of hashtable
    N_HT = 3  # number of hashtables

    ht = khmer.new_hashbits(K, HT_SIZE, N_HT)
    ht.consume_fasta_and_tag(filename)

    subset = ht.do_subset_partition(0, 0)
    ht.merge_subset(subset)

    n, _ = ht.count_partitions()
    assert n == 49

    ht.find_unpart(filename2, False, False)     # <-- don't traverse
    n, _ = ht.count_partitions()
    assert n == 99, n                    # all sequences disconnected


@attr('highmem')
def test_find_unpart_fail():
    filename = utils.get_test_data('random-20-a.odd.fa')
    filename2 = utils.get_test_data('random-20-a.odd.fa')  # <- switch to odd

    K = 20  # size of kmer
    HT_SIZE = 1e4  # size of hashtable
    N_HT = 3  # number of hashtables

    ht = khmer.new_hashbits(K, HT_SIZE, N_HT)
    ht.consume_fasta_and_tag(filename)

    subset = ht.do_subset_partition(0, 0)
    ht.merge_subset(subset)

    n, _ = ht.count_partitions()
    assert n == 49

    ht.find_unpart(filename2, True, False)
    n, _ = ht.count_partitions()
    assert n == 49, n                    # only 49 sequences worth of tags


def test_simple_median():
    hi = khmer.new_hashbits(6, 2, 2)

    (median, average, stddev) = hi.get_median_count("AAAAAA")
    print median, average, stddev
    assert median == 0
    assert average == 0.0
    assert stddev == 0.0

    hi.consume("AAAAAA")
    (median, average, stddev) = hi.get_median_count("AAAAAA")
    print median, average, stddev
    assert median == 1
    assert average == 1.0
    assert stddev == 0.0

####


def test_load_notexist_should_fail():
    savepath = utils.get_temp_filename('temphashbitssave0.ht')

    hi = khmer.new_counting_hash(12, 2)
    try:
        hi.load(savepath)
        assert 0, "load should fail"
    except IOError:
        pass


def test_load_truncated_should_fail():
    inpath = utils.get_test_data('random-20-a.fa')
    savepath = utils.get_temp_filename('temphashbitssave0.kh')

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
    except IOError, e:
        print str(e)


def test_save_load_tagset_notexist():
    ht = khmer.new_hashbits(32, 1, 1)

    outfile = utils.get_temp_filename('tagset')
    try:
        ht.load_tagset(outfile)
        assert 0, "this test should fail"
    except IOError, e:
        print str(e)


def test_save_load_tagset_trunc():
    ht = khmer.new_hashbits(32, 1, 1)

    outfile = utils.get_temp_filename('tagset')

    ht.add_tag('A' * 32)
    ht.add_tag('G' * 32)
    ht.save_tagset(outfile)
    ht.save_tagset('/tmp/goodversion-k32.tagset')

    # truncate tagset file...
    fp = open(outfile, 'rb')
    data = fp.read()
    fp.close()

    fp = open(outfile, 'wb')
    fp.write(data[:26])
    fp.close()

    # try loading it...
    try:
        ht.load_tagset(outfile)
        assert 0, "this test should fail"
    except IOError:
        pass

# to build the test files used below, add 'test' to this function
# and then look in /tmp. You will need to tweak the version info in
# khmer.hh in order to create "bad" versions, of course. -CTB


def _build_testfiles():
    # hashbits file

    inpath = utils.get_test_data('random-20-a.fa')
    hi = khmer.new_hashbits(12, 2)
    hi.consume_fasta(inpath)
    hi.save('/tmp/goodversion-k12.ht')

    # tagset file

    ht = khmer.new_hashbits(32, 1, 1)

    ht.add_tag('A' * 32)
    ht.add_tag('G' * 32)
    ht.save_tagset('/tmp/goodversion-k32.tagset')

    # stoptags file

    fakelump_fa = utils.get_test_data('fakelump.fa')

    ht = khmer.new_hashbits(32, 4, 4)
    ht.consume_fasta_and_tag(fakelump_fa)

    subset = ht.do_subset_partition(0, 0)
    ht.merge_subset(subset)

    EXCURSION_DISTANCE = 40
    EXCURSION_KMER_THRESHOLD = 82
    EXCURSION_KMER_COUNT_THRESHOLD = 1
    counting = khmer.new_counting_hash(32, 4, 4)

    ht.repartition_largest_partition(None, counting,
                                     EXCURSION_DISTANCE,
                                     EXCURSION_KMER_THRESHOLD,
                                     EXCURSION_KMER_COUNT_THRESHOLD)

    ht.save_stop_tags('/tmp/goodversion-k32.stoptags')


def test_hashbits_file_version_check():
    ht = khmer.new_hashbits(12, 1, 1)

    inpath = utils.get_test_data('badversion-k12.ht')

    try:
        ht.load(inpath)
        assert 0, "this should fail"
    except IOError, e:
        print str(e)


def test_hashbits_file_type_check():
    kh = khmer.new_counting_hash(12, 1, 1)
    savepath = utils.get_temp_filename('tempcountingsave0.kh')
    kh.save(savepath)

    ht = khmer.new_hashbits(12, 1, 1)

    try:
        ht.load(savepath)
        assert 0, "this should fail"
    except IOError, e:
        print str(e)


def test_stoptags_file_version_check():
    ht = khmer.new_hashbits(32, 1, 1)

    inpath = utils.get_test_data('badversion-k32.stoptags')

    try:
        ht.load_stop_tags(inpath)
        assert 0, "this should fail"
    except IOError, e:
        print str(e)


def test_stoptags_ksize_check():
    ht = khmer.new_hashbits(31, 1, 1)

    inpath = utils.get_test_data('goodversion-k32.stoptags')
    try:
        ht.load_stop_tags(inpath)
        assert 0, "this should fail"
    except IOError, e:
        print str(e)


def test_stop_tags_filetype_check():
    ht = khmer.new_hashbits(31, 1, 1)

    inpath = utils.get_test_data('goodversion-k32.tagset')
    try:
        ht.load_stop_tags(inpath)
        assert 0, "this should fail"
    except IOError, e:
        print str(e)


def test_tagset_file_version_check():
    ht = khmer.new_hashbits(32, 1, 1)

    inpath = utils.get_test_data('badversion-k32.tagset')

    try:
        ht.load_tagset(inpath)
        assert 0, "this should fail"
    except IOError, e:
        print str(e)


def test_tagset_ksize_check():
    ht = khmer.new_hashbits(31, 1, 1)

    inpath = utils.get_test_data('goodversion-k32.tagset')
    try:
        ht.load_tagset(inpath)
        assert 0, "this should fail"
    except IOError, e:
        print str(e)


def test_tagset_filetype_check():
    ht = khmer.new_hashbits(31, 1, 1)

    inpath = utils.get_test_data('goodversion-k32.stoptags')
    try:
        ht.load_tagset(inpath)
        assert 0, "this should fail"
    except IOError, e:
        print str(e)


def test_bad_primes_list():
    try:
        coutingtable = khmer._new_hashbits(31, ["a", "b", "c"], 1)
        assert 0, "Bad primes list should fail"
    except TypeError, e:
        print str(e)


def test_consume_absentfasta_with_reads_parser():
    presencetable = khmer.new_hashbits(31, 1, 1)
    try:
        presencetable.consume_fasta_with_reads_parser()
        assert 0, "this should fail"
    except TypeError, err:
        print str(err)
    readparser = ReadParser(utils.get_test_data('empty-file'))
    try:
        presencetable.consume_fasta_with_reads_parser(readparser)
        assert 0, "this should fail"
    except IOError, err:
        print str(err)
