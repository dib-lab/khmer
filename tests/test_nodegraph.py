from __future__ import print_function
from __future__ import absolute_import
#
# This file is part of khmer, htabletps://github.com/dib-lab/khmer/, and is
# Copyrightable (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# pylint: disable=missing-docstring,protected-access,no-member,
import khmer
from khmer import ReadParser

from screed.fasta import fasta_iter
import screed

from . import khmer_tst_utils as utils
from nose.plugins.attrib import attr


def teardown():
    utils.cleanup()


@attr('huge')
def test_toobig():
    try:
        pt = khmer.Nodegraph(32, 1e13, 1)
        assert 0, "This should fail"
    except MemoryError as err:
        print(str(err))


def test__get_set_tag_density():
    htableable = khmer._Nodegraph(32, [1])

    orig = htableable._get_tag_density()
    assert orig != 2
    htableable._set_tag_density(2)
    assert htableable._get_tag_density() == 2


def test_update_from():
    htableable = khmer.Nodegraph(5, 1000, 4)
    other_htableable = khmer.Nodegraph(5, 1000, 4)

    assert htableable.get('AAAAA') == 0
    assert htableable.get('GCGCG') == 0
    assert other_htableable.get('AAAAA') == 0
    assert other_htableable.get('GCGCG') == 0

    other_htableable.count('AAAAA')

    assert htableable.get('AAAAA') == 0
    assert htableable.get('GCGCG') == 0
    assert other_htableable.get('AAAAA') == 1
    assert other_htableable.get('GCGCG') == 0

    htableable.count('GCGCG')

    assert htableable.get('AAAAA') == 0
    assert htableable.get('GCGCG') == 1
    assert other_htableable.get('AAAAA') == 1
    assert other_htableable.get('GCGCG') == 0

    htableable.update(other_htableable)

    assert htableable.get('AAAAA') == 1
    assert htableable.get('GCGCG') == 1
    assert other_htableable.get('AAAAA') == 1
    assert other_htableable.get('GCGCG') == 0


def test_update_from_diff_ksize_2():
    htableable = khmer.Nodegraph(5, 1000, 4)
    other_htableable = khmer.Nodegraph(4, 1000, 4)

    try:
        htableable.update(other_htableable)
        assert 0, "should not be reached"
    except ValueError as err:
        print(str(err))

    try:
        other_htableable.update(htableable)
        assert 0, "should not be reached"
    except ValueError as err:
        print(str(err))


def test_update_from_diff_tablesize():
    htableable = khmer.Nodegraph(5, 100, 4)
    other_htableable = khmer.Nodegraph(5, 1000, 4)

    try:
        htableable.update(other_htableable)
        assert 0, "should not be reached"
    except ValueError as err:
        print(str(err))


def test_update_from_diff_num_tables():
    htableable = khmer.Nodegraph(5, 1000, 3)
    other_htableable = khmer.Nodegraph(5, 1000, 4)

    try:
        htableable.update(other_htableable)
        assert 0, "should not be reached"
    except ValueError as err:
        print(str(err))


def test_n_occupied_1():
    filename = utils.get_test_data('random-20-a.fa')

    ksize = 20  # size of kmer
    htable_size = 100000  # size of hashtableable
    num_htableables = 1  # number of hashtableables

    # test modified c++ n_occupied code
    htableable = khmer.Nodegraph(ksize, htable_size, num_htableables)

    for _, record in enumerate(fasta_iter(open(filename))):
        htableable.consume(record['sequence'])

    # this number calculated independently
    assert htableable.n_occupied() == 3884, htableable.n_occupied()


def test_bloom_python_1():
    # test python code to count unique kmers using bloom filter
    filename = utils.get_test_data('random-20-a.fa')

    ksize = 20  # size of kmer
    htable_size = 100000  # size of hashtableable
    num_htableables = 3  # number of hashtableables

    htableable = khmer.Nodegraph(ksize, htable_size, num_htableables)

    n_unique = 0
    for _, record in enumerate(fasta_iter(open(filename))):
        sequence = record['sequence']
        seq_len = len(sequence)
        for n in range(0, seq_len + 1 - ksize):
            kmer = sequence[n:n + ksize]
            if not htableable.get(kmer):
                n_unique += 1
            htableable.count(kmer)

    assert n_unique == 3960
    assert htableable.n_occupied() == 3885, htableable.n_occupied()

    # this number equals n_unique
    assert htableable.n_unique_kmers() == 3960, htableable.n_unique_kmers()


def test_bloom_c_1():
    # test c++ code to count unique kmers using bloom filter

    filename = utils.get_test_data('random-20-a.fa')

    ksize = 20  # size of kmer
    htable_size = 100000  # size of hashtableable
    num_htableables = 3  # number of hashtableables

    htableable = khmer.Nodegraph(ksize, htable_size, num_htableables)

    for _, record in enumerate(fasta_iter(open(filename))):
        htableable.consume(record['sequence'])

    assert htableable.n_occupied() == 3885
    assert htableable.n_unique_kmers() == 3960


def test_n_occupied_2():  # simple one
    ksize = 4
    htable_size = 10  # use 11
    num_htableables = 1

    htableable = khmer._Nodegraph(ksize, [11])
    htableable.count('AAAA')  # 00 00 00 00 = 0
    assert htableable.n_occupied() == 1

    htableable.count('ACTG')  # 00 10 01 11 =
    assert htableable.n_occupied() == 2

    htableable.count('AACG')  # 00 00 10 11 = 11  # collision 1

    assert htableable.n_occupied() == 2
    htableable.count('AGAC')   # 00  11 00 10 # collision 2
    assert htableable.n_occupied() == 2, htableable.n_occupied()


def test_bloom_c_2():  # simple one
    ksize = 4

    # use only 1 hashtableable, no bloom filter
    htableable = khmer._Nodegraph(ksize, [11])
    htableable.count('AAAA')  # 00 00 00 00 = 0
    htableable.count('ACTG')  # 00 10 01 11 =
    assert htableable.n_unique_kmers() == 2
    htableable.count('AACG')  # 00 00 10 11 = 11  # collision  with 1st kmer
    assert htableable.n_unique_kmers() == 2
    htableable.count('AGAC')   # 00  11 00 10 # collision  with 2nd kmer
    assert htableable.n_unique_kmers() == 2

    # use two hashtableables with 11,13
    other_htableable = khmer._Nodegraph(ksize, [11, 13])
    other_htableable.count('AAAA')  # 00 00 00 00 = 0

    other_htableable.count('ACTG')  # 00 10 01 11 = 2*16 +4 +3 = 39
    assert other_htableable.n_unique_kmers() == 2
    # 00 00 10 11 = 11  # collision with only 1st kmer
    other_htableable.count('AACG')
    assert other_htableable.n_unique_kmers() == 3
    other_htableable.count('AGAC')
    # 00  11 00 10  3*16 +2 = 50
    # collision with both 2nd and 3rd kmers

    assert other_htableable.n_unique_kmers() == 3


def test_filter_if_present():
    htable = khmer._Nodegraph(32, [3, 5])

    maskfile = utils.get_test_data('filter-test-A.fa')
    inputfile = utils.get_test_data('filter-test-B.fa')
    outfile = utils.get_temp_filename('filter')

    htable.consume_fasta(maskfile)
    htable.filter_if_present(inputfile, outfile)

    records = list(fasta_iter(open(outfile)))
    assert len(records) == 1
    assert records[0]['name'] == '3'


def test_combine_pe():
    inpfile = utils.get_test_data('combine_parts_1.fa')
    htable = khmer._Nodegraph(32, [1])

    htable.consume_partitioned_fasta(inpfile)
    assert htable.count_partitions() == (2, 0)

    first_seq = "CATGCAGAAGTTCCGCAACCATACCGTTCAGT"
    pid1 = htable.get_partition_id(first_seq)

    second_seq = "CAAATGTACATGCACTTAAAATCATCCAGCCG"
    pid2 = htable.get_partition_id(second_seq)

    assert pid1 == 2
    assert pid2 == 80293

    htable.join_partitions(pid1, pid2)

    pid1 = htable.get_partition_id(first_seq)
    pid2 = htable.get_partition_id(second_seq)

    assert pid1 == pid2
    assert htable.count_partitions() == (1, 0)


def test_load_partitioned():
    inpfile = utils.get_test_data('combine_parts_1.fa')
    htable = khmer._Nodegraph(32, [1])

    htable.consume_partitioned_fasta(inpfile)
    assert htable.count_partitions() == (2, 0)

    first_seq = "CATGCAGAAGTTCCGCAACCATACCGTTCAGT"
    assert htable.get(first_seq)

    second_seq = "CAAATGTACATGCACTTAAAATCATCCAGCCG"
    assert htable.get(second_seq)

    third_s = "CATGCAGAAGTTCCGCAACCATACCGTTCAGTTCCTGGTGGCTA"[-32:]
    assert htable.get(third_s)


def test_count_within_radius_simple():
    inpfile = utils.get_test_data('all-A.fa')
    htable = khmer._Nodegraph(4, [3, 5])

    print(htable.consume_fasta(inpfile))
    n = htable.count_kmers_within_radius('AAAA', 1)
    assert n == 1

    n = htable.count_kmers_within_radius('AAAA', 10)
    assert n == 1


def test_count_within_radius_big():
    inpfile = utils.get_test_data('random-20-a.fa')
    htable = khmer.Nodegraph(20, 1e5, 4)

    htable.consume_fasta(inpfile)
    n = htable.count_kmers_within_radius('CGCAGGCTGGATTCTAGAGG', int(1e6))
    assert n == 3961, n

    htable = khmer.Nodegraph(21, 1e5, 4)
    htable.consume_fasta(inpfile)
    n = htable.count_kmers_within_radius('CGCAGGCTGGATTCTAGAGGC', int(1e6))
    assert n == 39


def test_count_kmer_degree():
    inpfile = utils.get_test_data('all-A.fa')
    htable = khmer._Nodegraph(4, [3, 5])
    htable.consume_fasta(inpfile)

    assert htable.kmer_degree('AAAA') == 2
    assert htable.kmer_degree('AAAT') == 1
    assert htable.kmer_degree('AATA') == 0
    assert htable.kmer_degree('TAAA') == 1


def test_save_load_tagset():
    htable = khmer._Nodegraph(32, [1])

    outfile = utils.get_temp_filename('tagset')

    htable.add_tag('A' * 32)
    htable.save_tagset(outfile)

    htable.add_tag('G' * 32)

    htable.load_tagset(outfile)              # implicitly => clear_tags=True
    htable.save_tagset(outfile)

    # if tags have been cleared, then the new tagfile will be larger (34 bytes)
    # else smaller (26 bytes).

    fp = open(outfile, 'rb')
    data = fp.read()
    fp.close()
    assert len(data) == 30, len(data)


def test_save_load_tagset_noclear():
    htable = khmer._Nodegraph(32, [1])

    outfile = utils.get_temp_filename('tagset')

    htable.add_tag('A' * 32)
    htable.save_tagset(outfile)

    htable.add_tag('G' * 32)

    htable.load_tagset(outfile, False)  # set clear_tags => False; zero tags
    htable.save_tagset(outfile)

    # if tags have been cleared, then the new tagfile will be large (34 bytes);
    # else small (26 bytes).

    fp = open(outfile, 'rb')
    data = fp.read()
    fp.close()
    assert len(data) == 38, len(data)


def test_stop_traverse():
    filename = utils.get_test_data('random-20-a.fa')

    ksize = 20  # size of kmer
    htable_size = 1e4  # size of hashtableable
    num_htableables = 3  # number of hashtableables

    htable = khmer.Nodegraph(ksize, htable_size, num_htableables)

    # without tagging/joining across consume, this breaks into two partition;
    # with, it is one partition.
    htable.add_stop_tag('TTGCATACGTTGAGCCAGCG')

    # DO NOT join reads across stoptags
    htable.consume_fasta_and_tag(filename)
    subset = htable.do_subset_partition(0, 0, True)
    htable.merge_subset(subset)

    n, _ = htable.count_partitions()
    assert n == 2, n


def test_tag_across_stoptraverse():
    filename = utils.get_test_data('random-20-a.fa')

    ksize = 20  # size of kmer
    htable_size = 1e4  # size of hashtableable
    num_htableables = 3  # number of hashtableables

    htable = khmer.Nodegraph(ksize, htable_size, num_htableables)

    # without tagging/joining across consume, this breaks into two partition;
    # with, it is one partition.
    htable.add_stop_tag('CCGAATATATAACAGCGACG')

    # DO join reads across
    htable.consume_fasta_and_tag_with_stoptags(filename)
    subset = htable.do_subset_partition(0, 0)
    n, _ = htable.count_partitions()
    assert n == 99                       # reads only connected by traversal...

    n, _ = htable.subset_count_partitions(subset)
    assert n == 2                        # but need main to cross stoptags.

    htable.merge_subset(subset)

    n, _ = htable.count_partitions()         # ta-da!
    assert n == 1, n


def test_notag_across_stoptraverse():
    filename = utils.get_test_data('random-20-a.fa')

    ksize = 20  # size of kmer
    htable_size = 1e4  # size of hashtableable
    num_htableables = 3  # number of hashtableables

    htable = khmer.Nodegraph(ksize, htable_size, num_htableables)

    # connecting k-mer at the beginning/end of a read: breaks up into two.
    htable.add_stop_tag('TTGCATACGTTGAGCCAGCG')

    htable.consume_fasta_and_tag_with_stoptags(filename)

    subset = htable.do_subset_partition(0, 0)
    htable.merge_subset(subset)

    n, _ = htable.count_partitions()
    assert n == 2, n


def test_find_stoptags():
    htable = khmer._Nodegraph(5, [1])
    htable.add_stop_tag("AAAAA")

    assert htable.identify_stoptags_by_position("AAAAA") == [0]
    assert htable.identify_stoptags_by_position("AAAAAA") == [0, 1]
    assert htable.identify_stoptags_by_position("TTTTT") == [0]
    assert htable.identify_stoptags_by_position("TTTTTT") == [0, 1]


def test_find_stoptagsecond_seq():
    htable = khmer._Nodegraph(4, [1])
    htable.add_stop_tag("ATGC")

    x = htable.identify_stoptags_by_position("ATGCATGCGCAT")
    assert x == [0, 2, 4, 8], x


def test_get_ksize():
    kh = khmer._Nodegraph(22, [1])
    assert kh.ksize() == 22


def test_get_hashsizes():
    kh = khmer.Nodegraph(22, 100, 4)
    # Py2/3 hack, longify converts to long in py2, remove once py2 isn't
    # supported any longer.
    expected = utils.longify([97, 89, 83, 79])
    assert kh.hashsizes() == expected, kh.hashsizes()


def test_extract_unique_paths_0():
    kh = khmer._Nodegraph(10, [5, 7, 11, 13])

    x = kh.extract_unique_paths('ATGGAGAGACACAGATAGACAGGAGTGGCGATG', 10, 1)
    assert x == ['ATGGAGAGACACAGATAGACAGGAGTGGCGATG']

    kh.consume('ATGGAGAGACACAGATAGACAGGAGTGGCGATG')
    x = kh.extract_unique_paths('ATGGAGAGACACAGATAGACAGGAGTGGCGATG', 10, 1)
    assert not x


def test_extract_unique_paths_1():
    kh = khmer._Nodegraph(10, [5, 7, 11, 13])

    kh.consume('AGTGGCGATG')
    x = kh.extract_unique_paths('ATGGAGAGACACAGATAGACAGGAGTGGCGATG', 10, 1)
    print(x)
    assert x == ['ATGGAGAGACACAGATAGACAGGAGTGGCGAT']  # all but the last k-mer


def test_extract_unique_paths_2():
    kh = khmer._Nodegraph(10, [5, 7, 11, 13])

    kh.consume('ATGGAGAGAC')
    x = kh.extract_unique_paths('ATGGAGAGACACAGATAGACAGGAGTGGCGATG', 10, 1)
    print(x)
    assert x == ['TGGAGAGACACAGATAGACAGGAGTGGCGATG']  # all but the 1st k-mer


def test_extract_unique_paths_3():
    kh = khmer._Nodegraph(10, [5, 7, 11, 13])

    kh.consume('ATGGAGAGAC')
    kh.consume('AGTGGCGATG')
    x = kh.extract_unique_paths('ATGGAGAGACACAGATAGACAGGAGTGGCGATG', 10, 1)
    print(x)
    # all but the 1st/last k-mer
    assert x == ['TGGAGAGACACAGATAGACAGGAGTGGCGAT']


def test_extract_unique_paths_4():
    kh = khmer.Nodegraph(10, 1e6, 4)

    kh.consume('ATGGAGAGAC')
    kh.consume('AGTGGCGATG')

    kh.consume('ATAGACAGGA')

    x = kh.extract_unique_paths('ATGGAGAGACACAGATAGACAGGAGTGGCGATG', 10, 1)
    print(x)
    assert x == ['TGGAGAGACACAGATAGACAGG', 'TAGACAGGAGTGGCGAT']


def test_find_unpart():
    filename = utils.get_test_data('random-20-a.odd.fa')
    filename2 = utils.get_test_data('random-20-a.even.fa')

    ksize = 20  # size of kmer
    htable_size = 1e4  # size of hashtableable
    num_htableables = 3  # number of hashtableables

    htable = khmer.Nodegraph(ksize, htable_size, num_htableables)
    htable.consume_fasta_and_tag(filename)

    subset = htable.do_subset_partition(0, 0)
    htable.merge_subset(subset)

    n, _ = htable.count_partitions()
    assert n == 49

    htable.find_unpart(filename2, True, False)
    n, _ = htable.count_partitions()
    assert n == 1, n                     # all sequences connect


def test_find_unpart_notraverse():
    filename = utils.get_test_data('random-20-a.odd.fa')
    filename2 = utils.get_test_data('random-20-a.even.fa')

    ksize = 20  # size of kmer
    htable_size = 1e4  # size of hashtableable
    num_htableables = 3  # number of hashtableables

    htable = khmer.Nodegraph(ksize, htable_size, num_htableables)
    htable.consume_fasta_and_tag(filename)

    subset = htable.do_subset_partition(0, 0)
    htable.merge_subset(subset)

    n, _ = htable.count_partitions()
    assert n == 49

    htable.find_unpart(filename2, False, False)     # <-- don't traverse
    n, _ = htable.count_partitions()
    assert n == 99, n                    # all sequences disconnected


def test_find_unpart_fail():
    filename = utils.get_test_data('random-20-a.odd.fa')
    filename2 = utils.get_test_data('random-20-a.odd.fa')  # <- switch to odd

    ksize = 20  # size of kmer
    htable_size = 1e4  # size of hashtableable
    num_htableables = 3  # number of hashtableables

    htable = khmer.Nodegraph(ksize, htable_size, num_htableables)
    htable.consume_fasta_and_tag(filename)

    subset = htable.do_subset_partition(0, 0)
    htable.merge_subset(subset)

    n, _ = htable.count_partitions()
    assert n == 49

    htable.find_unpart(filename2, True, False)
    n, _ = htable.count_partitions()
    assert n == 49, n                    # only 49 sequences worth of tags


def test_simple_median():
    hi = khmer.Nodegraph(6, 1e5, 2)

    (median, average, stddev) = hi.get_median_count("AAAAAA")
    print(median, average, stddev)
    assert median == 0
    assert average == 0.0
    assert stddev == 0.0

    hi.consume("AAAAAA")
    (median, average, stddev) = hi.get_median_count("AAAAAA")
    print(median, average, stddev)
    assert median == 1
    assert average == 1.0
    assert stddev == 0.0


def test_badget():
    hbts = khmer.Nodegraph(6, 1e6, 1)

    dna = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAG"

    hbts.consume(dna)

    assert hbts.get("AGCTTT") == 1

    assert hbts.get("GATGAG") == 0

    try:
        hbts.get(b"AGCTT")
        assert 0, "this should fail"
    except ValueError as err:
        print(str(err))

    try:
        hbts.get(u"AGCTT")
        assert 0, "this should fail"
    except ValueError as err:
        print(str(err))


#


def test_load_notexist_should_fail():
    savepath = utils.get_temp_filename('tempnodegraphsave0.htable')

    hi = khmer._Countgraph(12, [1])
    try:
        hi.load(savepath)
        assert 0, "load should fail"
    except OSError:
        pass


def test_load_truncated_should_fail():
    inpath = utils.get_test_data('random-20-a.fa')
    savepath = utils.get_temp_filename('tempnodegraphsave0.ct')

    hi = khmer.Countgraph(12, 1000, 2)

    hi.consume_fasta(inpath)
    hi.save(savepath)

    fp = open(savepath, 'rb')
    data = fp.read()
    fp.close()

    fp = open(savepath, 'wb')
    fp.write(data[:1000])
    fp.close()

    hi = khmer._Countgraph(12, [1])
    try:
        hi.load(savepath)
        assert 0, "load should fail"
    except OSError as e:
        print(str(e))


def test_save_load_tagset_notexist():
    htable = khmer._Nodegraph(32, [1])

    outfile = utils.get_temp_filename('tagset')
    try:
        htable.load_tagset(outfile)
        assert 0, "this test should fail"
    except OSError as e:
        print(str(e))


def test_save_load_tagset_trunc():
    htable = khmer._Nodegraph(32, [1])

    outfile = utils.get_temp_filename('tagset')

    htable.add_tag('A' * 32)
    htable.add_tag('G' * 32)
    htable.save_tagset(outfile)

    # truncate tagset file...
    fp = open(outfile, 'rb')
    data = fp.read()
    fp.close()

    for i in range(len(data)):
        fp = open(outfile, 'wb')
        fp.write(data[:i])
        fp.close()

        # try loading it...
        try:
            htable.load_tagset(outfile)
            assert 0, "this test should fail"
        except OSError as err:
            print(str(err), i)

    # try loading it...
    try:
        htable.load_tagset(outfile)
        assert 0, "this test should fail"
    except OSError:
        pass

# to build the test files used below, add 'test' to this function
# and then look in /tmp. You will need to tweak the version info in
# khmer.hh in order to create "bad" versions, of course. -CTB


def _build_testfiles():
    # nodegraph file

    inpath = utils.get_test_data('random-20-a.fa')
    hi = khmer.Nodegraph(12, 2)
    hi.consume_fasta(inpath)
    hi.save('/tmp/goodversion-k12.htable')

    # tagset file

    htable = khmer._Nodegraph(32, [1])

    htable.add_tag('A' * 32)
    htable.add_tag('G' * 32)
    htable.save_tagset('/tmp/goodversion-k32.tagset')

    # stoptags file

    fakelump_fa = utils.get_test_data('fakelump.fa')

    htable = khmer.Nodegraph(32, 4, 4)
    htable.consume_fasta_and_tag(fakelump_fa)

    subset = htable.do_subset_partition(0, 0)
    htable.merge_subset(subset)

    EXCURSION_DISTANCE = 40
    EXCURSION_ksizeMER_THRESHOLD = 82
    EXCURSION_ksizeMER_COUNT_THRESHOLD = 1
    counting = khmer.Countgraph(32, 4, 4)

    htable.repartition_largest_partition(None, counting,
                                         EXCURSION_DISTANCE,
                                         EXCURSION_ksizeMER_THRESHOLD,
                                         EXCURSION_ksizeMER_COUNT_THRESHOLD)

    htable.save_stop_tags('/tmp/goodversion-k32.stoptags')


def test_nodegraph_file_version_check():
    htable = khmer._Nodegraph(12, [1])

    inpath = utils.get_test_data('badversion-k12.htable')

    try:
        htable.load(inpath)
        assert 0, "this should fail"
    except OSError as e:
        print(str(e))


def test_nodegraph_file_type_check():
    kh = khmer._Countgraph(12, [1])
    savepath = utils.get_temp_filename('tempcountingsave0.ct')
    kh.save(savepath)

    htable = khmer._Nodegraph(12, [1])

    try:
        htable.load(savepath)
        assert 0, "this should fail"
    except OSError as e:
        print(str(e))


def test_stoptags_file_version_check():
    htable = khmer._Nodegraph(32, [1])

    inpath = utils.get_test_data('badversion-k32.stoptags')

    try:
        htable.load_stop_tags(inpath)
        assert 0, "this should fail"
    except OSError as e:
        print(str(e))


def test_stoptags_ksize_check():
    htable = khmer._Nodegraph(31, [1])

    inpath = utils.get_test_data('goodversion-k32.stoptags')
    try:
        htable.load_stop_tags(inpath)
        assert 0, "this should fail"
    except OSError as e:
        print(str(e))


def test_stop_tags_filetype_check():
    htable = khmer._Nodegraph(31, [1])

    inpath = utils.get_test_data('goodversion-k32.tagset')
    try:
        htable.load_stop_tags(inpath)
        assert 0, "this should fail"
    except OSError as e:
        print(str(e))


def test_tagset_file_version_check():
    htable = khmer._Nodegraph(32, [1])

    inpath = utils.get_test_data('badversion-k32.tagset')

    try:
        htable.load_tagset(inpath)
        assert 0, "this should fail"
    except OSError as e:
        print(str(e))


def test_stop_tags_truncate_check():
    htable = khmer._Nodegraph(32, [1])

    inpath = utils.get_test_data('goodversion-k32.tagset')
    data = open(inpath, 'rb').read()

    truncpath = utils.get_temp_filename('zzz')
    for i in range(len(data)):
        fp = open(truncpath, 'wb')
        fp.write(data[:i])
        fp.close()

        try:
            htable.load_stop_tags(truncpath)
            assert 0, "expect failure of previous command"
        except OSError as e:
            print(i, str(e))


def test_tagset_ksize_check():
    htable = khmer._Nodegraph(31, [1])

    inpath = utils.get_test_data('goodversion-k32.tagset')
    try:
        htable.load_tagset(inpath)
        assert 0, "this should fail"
    except OSError as e:
        print(str(e))


def test_tagset_filetype_check():
    htable = khmer._Nodegraph(31, [1])

    inpath = utils.get_test_data('goodversion-k32.stoptags')
    try:
        htable.load_tagset(inpath)
        assert 0, "this should fail"
    except OSError as e:
        print(str(e))


def test_bad_primes_list():
    try:
        coutingtable = khmer._Nodegraph(31, ["a", "b", "c"], 1)
        assert 0, "Bad primes list should fail"
    except TypeError as e:
        print(str(e))


def test_consume_absentfasta_with_reads_parser():
    presencetable = khmer._Nodegraph(31, [1])
    try:
        presencetable.consume_fasta_with_reads_parser()
        assert 0, "this should fail"
    except TypeError as err:
        print(str(err))
    try:
        readparser = ReadParser(utils.get_test_data('empty-file'))
        presencetable.consume_fasta_with_reads_parser(readparser)
        assert 0, "this should fail"
    except OSError as err:
        print(str(err))
    except ValueError as err:
        print(str(err))


def test_bad_primes():
    try:
        countingtable = khmer._Nodegraph.__new__(
            khmer._Nodegraph, 6, ["a", "b", "c"])
        assert 0, "this should fail"
    except TypeError as e:
        print(str(e))


def test_consume_fasta_and_tag_with_badreads_parser():
    presencetable = khmer.Nodegraph(6, 1e6, 2)
    try:
        readsparser = khmer.ReadParser(utils.get_test_data("test-empty.fa"))
        presencetable.consume_fasta_and_tag_with_reads_parser(readsparser)
        assert 0, "this should fail"
    except OSError as e:
        print(str(e))
    except ValueError as e:
        print(str(e))