# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2010-2015, Michigan State University.
# Copyright (C) 2015-2016, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org
# pylint: disable=missing-docstring,protected-access,no-member,invalid-name


import khmer
from khmer import Nodegraph, Countgraph
from khmer import ReadParser
from khmer import reverse_complement as revcomp
from khmer.khmer_args import create_matching_nodegraph

import screed

import pytest

from . import khmer_tst_utils as utils


def teardown():
    utils.cleanup()


@pytest.mark.huge
def test_toobig():
    try:
        khmer.Nodegraph(32, 1e13, 1)
        assert 0, "This should fail"
    except MemoryError as err:
        print(str(err))


def test_add_tag():
    nodegraph = khmer.Nodegraph(6, 1, 1)

    assert nodegraph.n_tags == 0
    nodegraph.add_tag('AATAAG')
    assert nodegraph.n_tags == 1

    print(nodegraph.get_tagset())
    assert nodegraph.get_tagset() == ['AATAAG']


def test_add_tag():
    nodegraph = khmer.Nodegraph(6, 1, 1)

    assert nodegraph.n_tags == 0
    nodegraph.add_tag('AATAAG')
    assert nodegraph.n_tags == 1

    print(nodegraph.get_tagset())
    assert nodegraph.get_tagset() == ['AATAAG']


def test_get_tag_sequences():
    nodegraph = khmer.Nodegraph(6, 1, 1)

    assert nodegraph.n_tags == 0
    kmer = nodegraph.hash('AATAAG')
    nodegraph.add_tag(kmer)
    assert nodegraph.n_tags == 1

    tagset = nodegraph.get_tags_for_sequence('GGGAATAAGGGG')
    tagset = list(tagset)

    assert len(tagset) == 1
    assert nodegraph.reverse_hash(tagset[0]) == 'AATAAG'


def test_add_stop_tag():
    nodegraph = khmer.Nodegraph(6, 1, 1)

    nodegraph.add_stop_tag('AATAAG')
    print(nodegraph.get_stop_tags())
    assert nodegraph.get_stop_tags() == ['AATAAG']


def test_add_stop_tag_hashval():
    nodegraph = khmer.Nodegraph(6, 1, 1)

    kmer = nodegraph.hash('AATAAG')
    nodegraph.add_stop_tag(kmer)
    print(nodegraph.get_stop_tags())
    assert nodegraph.get_stop_tags() == ['AATAAG']


def test__get_set_tag_density():
    nodegraph = khmer.Nodegraph(32, 1, 1,)
    orig = nodegraph.tag_density
    assert orig != 2
    nodegraph.tag_density = 2
    assert nodegraph.tag_density == 2


def test_update_from():
    nodegraph = khmer.Nodegraph(5, 1000, 4)
    other_nodegraph = khmer.Nodegraph(5, 1000, 4)

    assert nodegraph.get('AAAAA') == 0
    assert nodegraph.get('GCGCG') == 0
    assert nodegraph.n_occupied() == 0
    assert other_nodegraph.get('AAAAA') == 0
    assert other_nodegraph.get('GCGCG') == 0
    assert other_nodegraph.n_occupied() == 0

    other_nodegraph.count('AAAAA')

    assert nodegraph.get('AAAAA') == 0
    assert nodegraph.get('GCGCG') == 0
    assert nodegraph.n_occupied() == 0
    assert other_nodegraph.get('AAAAA') == 1
    assert other_nodegraph.get('GCGCG') == 0
    assert other_nodegraph.n_occupied() == 1

    nodegraph.count('GCGCG')

    assert nodegraph.get('AAAAA') == 0
    assert nodegraph.get('GCGCG') == 1
    assert nodegraph.n_occupied() == 1
    assert other_nodegraph.get('AAAAA') == 1
    assert other_nodegraph.get('GCGCG') == 0
    assert other_nodegraph.n_occupied() == 1

    nodegraph.update(other_nodegraph)

    assert nodegraph.get('AAAAA') == 1
    assert nodegraph.get('GCGCG') == 1
    assert nodegraph.n_occupied() == 2
    assert other_nodegraph.get('AAAAA') == 1
    assert other_nodegraph.get('GCGCG') == 0
    assert other_nodegraph.n_occupied() == 1


def test_update_from_2():

    ng1 = khmer.Nodegraph(20, 1000, 4)
    ng2 = khmer.Nodegraph(20, 1000, 4)

    filename = utils.get_test_data('random-20-a.fa')
    ng1.consume_seqfile(filename)
    ng2.consume_seqfile(filename)

    assert ng1.n_occupied() == ng2.n_occupied()
    ng1.update(ng2)

    assert ng1.n_occupied() == ng2.n_occupied()


def test_update_from_diff_ksize_2():
    nodegraph = khmer.Nodegraph(5, 1000, 4)
    other_nodegraph = khmer.Nodegraph(4, 1000, 4)

    try:
        nodegraph.update(other_nodegraph)
        assert 0, "should not be reached"
    except ValueError as err:
        print(str(err))

    try:
        other_nodegraph.update(nodegraph)
        assert 0, "should not be reached"
    except ValueError as err:
        print(str(err))


def test_update_from_diff_tablesize():
    nodegraph = khmer.Nodegraph(5, 100, 4)
    other_nodegraph = khmer.Nodegraph(5, 1000, 4)

    try:
        nodegraph.update(other_nodegraph)
        assert 0, "should not be reached"
    except ValueError as err:
        print(str(err))


def test_update_from_diff_num_tables():
    nodegraph = khmer.Nodegraph(5, 1000, 3)
    other_nodegraph = khmer.Nodegraph(5, 1000, 4)

    try:
        nodegraph.update(other_nodegraph)
        assert 0, "should not be reached"
    except ValueError as err:
        print(str(err))


def test_similarity_1():
    nodegraph = khmer.Nodegraph(5, 1000, 4)
    other_nodegraph = khmer.Nodegraph(5, 1000, 4)

    assert nodegraph.similarity(other_nodegraph) == 0

    other_nodegraph.count('AAAAA')

    assert nodegraph.similarity(other_nodegraph) == 0

    nodegraph.count('GCGCG')

    assert nodegraph.similarity(other_nodegraph) == 0

    nodegraph.count('AAAAA')
    other_nodegraph.count('GCGCG')

    assert nodegraph.similarity(other_nodegraph) == 1


def test_n_occupied_1():
    filename = utils.get_test_data('random-20-a.fa')

    ksize = 20  # size of kmer
    htable_size = 100000  # size of hashtable
    num_nodegraphs = 1  # number of hashtables

    # test modified c++ n_occupied code
    nodegraph = khmer.Nodegraph(ksize, htable_size, num_nodegraphs)

    for _, record in enumerate(screed.open(filename)):
        nodegraph.consume(record.sequence)

    # this number calculated independently
    assert nodegraph.n_occupied() == 3884, nodegraph.n_occupied()


def test_bloom_python_1():
    # test python code to count unique kmers using bloom filter
    filename = utils.get_test_data('random-20-a.fa')

    ksize = 20  # size of kmer
    htable_size = 100000  # size of hashtable
    num_nodegraphs = 3  # number of hashtables

    nodegraph = khmer.Nodegraph(ksize, htable_size, num_nodegraphs)

    n_unique = 0
    for _, record in enumerate(screed.open(filename)):
        sequence = record.sequence
        seq_len = len(sequence)
        for num in range(0, seq_len + 1 - ksize):
            kmer = sequence[num:num + ksize]
            if not nodegraph.get(kmer):
                n_unique += 1
            nodegraph.count(kmer)

    assert n_unique == 3960
    assert nodegraph.n_occupied() == 3884, nodegraph.n_occupied()

    # this number equals n_unique
    assert nodegraph.n_unique_kmers() == 3960, nodegraph.n_unique_kmers()


def test_bloom_c_1():
    # test c++ code to count unique kmers using bloom filter

    filename = utils.get_test_data('random-20-a.fa')

    ksize = 20  # size of kmer
    htable_size = 100000  # size of hashtable
    num_nodegraphs = 3  # number of hashtables

    nodegraph = khmer.Nodegraph(ksize, htable_size, num_nodegraphs)

    for _, record in enumerate(screed.open(filename)):
        nodegraph.consume(record.sequence)

    assert nodegraph.n_occupied() == 3884
    assert nodegraph.n_unique_kmers() == 3960


def test_n_occupied_2():  # simple one
    ksize = 4

    nodegraph = khmer.Nodegraph(ksize, 1, 1, primes=[11])
    nodegraph.count('AAAA')  # 00 00 00 00 = 0
    assert nodegraph.n_occupied() == 1

    nodegraph.count('ACTG')  # 00 10 01 11 =
    assert nodegraph.n_occupied() == 2

    nodegraph.count('AACG')  # 00 00 10 11 = 11  # collision 1

    assert nodegraph.n_occupied() == 2
    nodegraph.count('AGAC')   # 00  11 00 10 # collision 2
    assert nodegraph.n_occupied() == 2, nodegraph.n_occupied()


def test_n_occupied_2_add_is_count():  # 'add' synonym for 'count'
    ksize = 4

    nodegraph = khmer.Nodegraph(ksize, 1, 1, primes=[11])
    nodegraph.add('AAAA')  # 00 00 00 00 = 0
    assert nodegraph.n_occupied() == 1

    nodegraph.add('ACTG')  # 00 10 01 11 =
    assert nodegraph.n_occupied() == 2

    nodegraph.add('AACG')  # 00 00 10 11 = 11  # collision 1

    assert nodegraph.n_occupied() == 2
    nodegraph.add('AGAC')   # 00  11 00 10 # collision 2
    assert nodegraph.n_occupied() == 2, nodegraph.n_occupied()


def test_bloom_c_2():  # simple one
    ksize = 4

    # use only 1 hashtable, no bloom filter
    nodegraph = khmer.Nodegraph(ksize, 1, 1, primes=[11])
    nodegraph.count('AAAA')  # 00 00 00 00 = 0
    nodegraph.count('ACTG')  # 00 10 01 11 =
    assert nodegraph.n_unique_kmers() == 2
    nodegraph.count('AACG')  # 00 00 10 11 = 11  # collision  with 1st kmer
    assert nodegraph.n_unique_kmers() == 2
    nodegraph.count('AGAC')   # 00  11 00 10 # collision  with 2nd kmer
    assert nodegraph.n_unique_kmers() == 2

    # use two hashtables with 11,13
    other_nodegraph = khmer.Nodegraph(ksize, 1, 1, primes=[11, 13])
    other_nodegraph.count('AAAA')  # 00 00 00 00 = 0

    other_nodegraph.count('ACTG')  # 00 10 01 11 = 2*16 +4 +3 = 39
    assert other_nodegraph.n_unique_kmers() == 2
    # 00 00 10 11 = 11  # collision with only 1st kmer
    other_nodegraph.count('AACG')
    assert other_nodegraph.n_unique_kmers() == 3
    other_nodegraph.count('AGAC')
    # 00  11 00 10  3*16 +2 = 50
    # collision with both 2nd and 3rd kmers

    assert other_nodegraph.n_unique_kmers() == 3


def test_combine_pe():
    inpfile = utils.get_test_data('combine_parts_1.fa')
    nodegraph = khmer.Nodegraph(32, 1, 1)

    nodegraph.consume_partitioned_fasta(inpfile)
    assert nodegraph.count_partitions() == (2, 0)

    first_seq = "CATGCAGAAGTTCCGCAACCATACCGTTCAGT"
    pid1 = nodegraph.get_partition_id(first_seq)

    second_seq = "CAAATGTACATGCACTTAAAATCATCCAGCCG"
    pid2 = nodegraph.get_partition_id(second_seq)

    assert pid1 == 2
    assert pid2 == 80293

    nodegraph.join_partitions(pid1, pid2)

    pid1 = nodegraph.get_partition_id(first_seq)
    pid2 = nodegraph.get_partition_id(second_seq)

    assert pid1 == pid2
    assert nodegraph.count_partitions() == (1, 0)


def test_load_partitioned():
    inpfile = utils.get_test_data('combine_parts_1.fa')
    nodegraph = khmer.Nodegraph(32, 1, 1)

    nodegraph.consume_partitioned_fasta(inpfile)
    assert nodegraph.count_partitions() == (2, 0)

    first_seq = "CATGCAGAAGTTCCGCAACCATACCGTTCAGT"
    assert nodegraph.get(first_seq)

    second_seq = "CAAATGTACATGCACTTAAAATCATCCAGCCG"
    assert nodegraph.get(second_seq)

    third_s = "CATGCAGAAGTTCCGCAACCATACCGTTCAGTTCCTGGTGGCTA"[-32:]
    assert nodegraph.get(third_s)


def test_consume_partitioned_fail():
    inpfile = utils.get_test_data('test-reads.fa')
    nodegraph = khmer.Nodegraph(32, 1, 1)

    with pytest.raises(ValueError):
        nodegraph.consume_partitioned_fasta(inpfile)


def test_count_within_radius_simple():
    inpfile = utils.get_test_data('all-A.fa')
    nodegraph = khmer.Nodegraph(4, 1, 1, primes=[3, 5])

    print(nodegraph.consume_seqfile(inpfile))
    n = nodegraph.count_kmers_within_radius('AAAA', 1)
    assert n == 1

    n = nodegraph.count_kmers_within_radius('AAAA', 10)
    assert n == 1


def test_count_within_radius_big():
    inpfile = utils.get_test_data('random-20-a.fa')
    nodegraph = khmer.Nodegraph(20, 1e5, 4)

    nodegraph.consume_seqfile(inpfile)
    n = nodegraph.count_kmers_within_radius('CGCAGGCTGGATTCTAGAGG', int(1e6))
    assert n == 3961, n

    nodegraph = khmer.Nodegraph(21, 1e5, 4)
    nodegraph.consume_seqfile(inpfile)
    n = nodegraph.count_kmers_within_radius('CGCAGGCTGGATTCTAGAGGC', int(1e6))
    assert n == 39


def test_count_kmer_degree():
    inpfile = utils.get_test_data('all-A.fa')
    nodegraph = khmer.Nodegraph(4, 1, 1, primes=[3, 5])
    nodegraph.consume_seqfile(inpfile)

    assert nodegraph.kmer_degree('AAAA') == 2
    assert nodegraph.kmer_degree('AAAT') == 1
    assert nodegraph.kmer_degree('AATA') == 0
    assert nodegraph.kmer_degree('TAAA') == 1


def test_kmer_neighbors():
    inpfile = utils.get_test_data('all-A.fa')
    nodegraph = khmer.Nodegraph(4, 100, 1)
    nodegraph.consume_seqfile(inpfile)

    def n_to_str(x):
        return [str(i) for i in x]

    h = nodegraph.hash('AAAA')
    print(type('AAAA'))
    assert n_to_str(nodegraph.neighbors(
        h)) == ['AAAA', 'AAAA']       # AAAA on both sides
    assert n_to_str(nodegraph.neighbors('AAAA')) == [
        'AAAA', 'AAAA']  # AAAA on both sides

    h = nodegraph.hash('AAAT')
    assert n_to_str(nodegraph.neighbors(h)) == ['AAAA']       # AAAA on 1 side
    assert n_to_str(nodegraph.neighbors('AAAT')) == ['AAAA']  # AAAA on 1 side

    h = nodegraph.hash('AATA')
    assert nodegraph.neighbors(h) == []           # no neighbors
    assert n_to_str(nodegraph.neighbors('AATA')) == []      # AAAA on one side

    h = nodegraph.hash('TAAA')
    assert n_to_str(nodegraph.neighbors(h)) == ['AAAA']       # AAAA on both
    assert n_to_str(nodegraph.neighbors('TAAA')) == ['AAAA']  # AAAA on both


def test_kmer_neighbors_wrong_ksize():
    inpfile = utils.get_test_data('all-A.fa')
    nodegraph = khmer.Nodegraph(4, 1, 1, primes=[3, 5])
    nodegraph.consume_seqfile(inpfile)

    try:
        nodegraph.neighbors('AAAAA')
        assert 0, "neighbors() should fail with too long string"
    except ValueError:
        pass

    try:
        nodegraph.neighbors(b'AAAAA')
        assert 0, "neighbors() should fail with too long string"
    except ValueError:
        pass

    try:
        nodegraph.neighbors({})
        assert 0, "neighbors() should fail with non hash/str arg"
    except TypeError:
        pass


def test_save_load_tagset():
    nodegraph = khmer.Nodegraph(32, 1, 1)

    outfile = utils.get_temp_filename('tagset')

    nodegraph.add_tag('A' * 32)
    nodegraph.save_tagset(outfile)

    nodegraph.add_tag('G' * 32)

    nodegraph.load_tagset(outfile)              # implicitly => clear_tags=True
    nodegraph.save_tagset(outfile)

    # if tags have been cleared, then the new tagfile will be larger (34 bytes)
    # else smaller (26 bytes).

    fp = open(outfile, 'rb')
    data = fp.read()
    fp.close()
    assert len(data) == 30, len(data)


def test_save_load_tagset_noclear():
    nodegraph = khmer.Nodegraph(32, 1, 1)

    outfile = utils.get_temp_filename('tagset')

    nodegraph.add_tag('A' * 32)
    nodegraph.save_tagset(outfile)

    nodegraph.add_tag('G' * 32)

    nodegraph.load_tagset(outfile, False)  # set clear_tags => False; zero tags
    nodegraph.save_tagset(outfile)

    # if tags have been cleared, then the new tagfile will be large (34 bytes);
    # else small (26 bytes).

    fp = open(outfile, 'rb')
    data = fp.read()
    fp.close()
    assert len(data) == 38, len(data)


def test_stop_traverse():
    filename = utils.get_test_data('random-20-a.fa')

    ksize = 20  # size of kmer
    htable_size = 1e4  # size of hashtable
    num_nodegraphs = 3  # number of hashtables

    nodegraph = khmer.Nodegraph(ksize, htable_size, num_nodegraphs)

    # without tagging/joining across consume, this breaks into two partition;
    # with, it is one partition.
    nodegraph.add_stop_tag('TTGCATACGTTGAGCCAGCG')

    # DO NOT join reads across stoptags
    nodegraph.consume_seqfile_and_tag(filename)
    subset = nodegraph.do_subset_partition(0, 0, True)
    nodegraph.merge_subset(subset)

    n, _ = nodegraph.count_partitions()
    assert n == 2, n


def test_get_ksize():
    kh = khmer.Nodegraph(22, 1, 1)
    assert kh.ksize() == 22


def test_get_hashsizes():
    kh = khmer.Nodegraph(22, 100, 4)
    # Py2/3 hack, longify converts to long in py2, remove once py2 isn't
    # supported any longer.
    expected = utils.longify([97, 89, 83, 79])
    assert kh.hashsizes() == expected, kh.hashsizes()


def test_extract_unique_paths_0():
    kh = khmer.Nodegraph(10, 1, 1, primes=[5, 7, 11, 13])

    x = kh.extract_unique_paths('ATGGAGAGACACAGATAGACAGGAGTGGCGATG', 10, 1)
    assert x == ['ATGGAGAGACACAGATAGACAGGAGTGGCGATG']

    kh.consume('ATGGAGAGACACAGATAGACAGGAGTGGCGATG')
    x = kh.extract_unique_paths('ATGGAGAGACACAGATAGACAGGAGTGGCGATG', 10, 1)
    assert not x


def test_extract_unique_paths_1():
    kh = khmer.Nodegraph(10, 1, 1, primes=[5, 7, 11, 13])

    kh.consume('AGTGGCGATG')
    x = kh.extract_unique_paths('ATGGAGAGACACAGATAGACAGGAGTGGCGATG', 10, 1)
    print(x)
    assert x == ['ATGGAGAGACACAGATAGACAGGAGTGGCGAT']  # all but the last k-mer


def test_extract_unique_paths_2():
    kh = khmer.Nodegraph(10, 1, 1, primes=[5, 7, 11, 13])

    kh.consume('ATGGAGAGAC')
    x = kh.extract_unique_paths('ATGGAGAGACACAGATAGACAGGAGTGGCGATG', 10, 1)
    print(x)
    assert x == ['TGGAGAGACACAGATAGACAGGAGTGGCGATG']  # all but the 1st k-mer


def test_extract_unique_paths_3():
    kh = khmer.Nodegraph(10, 1, 1, primes=[5, 7, 11, 13])

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


def test_get_raw_tables():
    kh = khmer.Nodegraph(10, 1e6, 4)
    kh.consume('ATGGAGAGAC')
    kh.consume('AGTGGCGATG')
    kh.consume('ATAGACAGGA')
    tables = kh.get_raw_tables()

    for size, table in zip(kh.hashsizes(), tables):
        assert isinstance(table, memoryview)
        assert size == len(table)


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

    try:
        hi = Countgraph.load(savepath)
        assert 0, "load should fail"
    except OSError:
        pass


def test_load_truncated_should_fail():
    inpath = utils.get_test_data('random-20-a.fa')
    savepath = utils.get_temp_filename('tempnodegraphsave0.ct')

    hi = khmer.Countgraph(12, 1000, 2)

    hi.consume_seqfile(inpath)
    hi.save(savepath)

    fp = open(savepath, 'rb')
    data = fp.read()
    fp.close()

    fp = open(savepath, 'wb')
    fp.write(data[:1000])
    fp.close()

    try:
        hi = Countgraph.load(savepath)
        assert 0, "load should fail"
    except OSError as e:
        print(str(e))


def test_save_load_tagset_notexist():
    nodegraph = khmer.Nodegraph(32, 1, 1)

    outfile = utils.get_temp_filename('tagset')
    try:
        nodegraph.load_tagset(outfile)
        assert 0, "this test should fail"
    except OSError as e:
        print(str(e))


def test_save_load_tagset_trunc():
    nodegraph = khmer.Nodegraph(32, 1, 1)

    outfile = utils.get_temp_filename('tagset')

    nodegraph.add_tag('A' * 32)
    nodegraph.add_tag('G' * 32)
    nodegraph.save_tagset(outfile)

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
            nodegraph.load_tagset(outfile)
            assert 0, "this test should fail"
        except OSError as err:
            print(str(err), i)

    # try loading it...
    try:
        nodegraph.load_tagset(outfile)
        assert 0, "this test should fail"
    except OSError:
        pass

# to build the test files used below, add 'test' to this function
# and then look in /tmp. You will need to tweak the version info in
# khmer.hh in order to create "bad" versions, of course. -CTB


def _build_testfiles():
    # nodegraph file

    inpath = utils.get_test_data('random-20-a.fa')
    hi = khmer.Nodegraph(12, 1, 1, primes=[2])
    hi.consume_seqfile(inpath)
    hi.save('/tmp/goodversion-k12.htable')

    # tagset file

    nodegraph = khmer.Nodegraph(32, 1, 1)

    nodegraph.add_tag('A' * 32)
    nodegraph.add_tag('G' * 32)
    nodegraph.save_tagset('/tmp/goodversion-k32.tagset')

    # stoptags file

    fakelump_fa = utils.get_test_data('fakelump.fa')

    nodegraph = khmer.Nodegraph(32, 100000, 4)
    nodegraph.consume_seqfile_and_tag(fakelump_fa)

    subset = nodegraph.do_subset_partition(0, 0)
    nodegraph.merge_subset(subset)

    EXCURSION_DISTANCE = 40
    EXCURSION_KMER_THRESHOLD = 82
    EXCURSION_KMER_COUNT_THRESHOLD = 1
    counting = khmer.Countgraph(32, 100000, 4)

    nodegraph.repartition_largest_partition(counting,
                                            EXCURSION_DISTANCE,
                                            EXCURSION_KMER_THRESHOLD,
                                            EXCURSION_KMER_COUNT_THRESHOLD)

    nodegraph.save_stop_tags('/tmp/goodversion-k32.stoptags')


def test_hashbits_file_version_check():

    inpath = utils.get_test_data('badversion-k12.htable')

    try:
        nodegraph = Nodegraph.load(inpath)
        assert 0, "this should fail"
    except OSError as e:
        print(str(e))


def test_nodegraph_file_type_check():
    kh = khmer.Countgraph(12, 1, 1)
    savepath = utils.get_temp_filename('tempcountingsave0.ct')
    kh.save(savepath)

    try:
        nodegraph = Nodegraph.load(savepath)
        assert 0, "this should fail"
    except OSError as e:
        print(str(e))


def test_stoptags_file_version_check():
    nodegraph = khmer.Nodegraph(32, 1, 1,)

    inpath = utils.get_test_data('badversion-k32.stoptags')

    try:
        nodegraph.load_stop_tags(inpath)
        assert 0, "this should fail"
    except OSError as e:
        print(str(e))


def test_stoptags_ksize_check():
    nodegraph = khmer.Nodegraph(31, 1, 1)

    inpath = utils.get_test_data('goodversion-k32.stoptags')
    try:
        nodegraph.load_stop_tags(inpath)
        assert 0, "this should fail"
    except OSError as e:
        print(str(e))


def test_stop_tags_filetype_check():
    nodegraph = khmer.Nodegraph(31, 1, 1)

    inpath = utils.get_test_data('goodversion-k32.tagset')
    try:
        nodegraph.load_stop_tags(inpath)
        assert 0, "this should fail"
    except OSError as e:
        print(str(e))


def test_tagset_file_version_check():
    nodegraph = khmer.Nodegraph(32, 1, 1)

    inpath = utils.get_test_data('badversion-k32.tagset')

    try:
        nodegraph.load_tagset(inpath)
        assert 0, "this should fail"
    except OSError as e:
        print(str(e))


def test_stop_tags_truncate_check():
    nodegraph = khmer.Nodegraph(32, 1, 1)

    inpath = utils.get_test_data('goodversion-k32.tagset')
    data = open(inpath, 'rb').read()

    truncpath = utils.get_temp_filename('zzz')
    for i in range(len(data)):
        fp = open(truncpath, 'wb')
        fp.write(data[:i])
        fp.close()

        try:
            nodegraph.load_stop_tags(truncpath)
            assert 0, "expect failure of previous command"
        except OSError as e:
            print(i, str(e))


def test_tagset_ksize_check():
    nodegraph = khmer.Nodegraph(31, 1, 1)

    inpath = utils.get_test_data('goodversion-k32.tagset')
    try:
        nodegraph.load_tagset(inpath)
        assert 0, "this should fail"
    except OSError as e:
        print(str(e))


def test_tagset_filetype_check():
    nodegraph = khmer.Nodegraph(31, 1, 1)

    inpath = utils.get_test_data('goodversion-k32.stoptags')
    try:
        nodegraph.load_tagset(inpath)
        assert 0, "this should fail"
    except OSError as e:
        print(str(e))


def test_bad_primes_list():
    try:
        khmer.Nodegraph(31, 1, 1, primes=["a", "b", "c"])
        assert 0, "Bad primes list should fail"
    except TypeError as e:
        print(str(e))


def test_consume_absentfasta():
    nodegraph = khmer.Nodegraph(31, 1, 1)
    try:
        nodegraph.consume_seqfile()
        assert 0, "this should fail"
    except TypeError as err:
        print(str(err))
    try:
        readparser = ReadParser(utils.get_test_data('empty-file'))
        nodegraph.consume_seqfile(readparser)
        assert 0, "this should fail"
    except OSError as err:
        print(str(err))
    except ValueError as err:
        print(str(err))


def test_bad_primes():
    try:
        Nodegraph(6, 1, 1, primes=["a", "b", "c"])
        assert 0, "this should fail"
    except TypeError as e:
        print(str(e))


def test_consume_seqfile_and_tag_with_badreads_parser():
    nodegraph = khmer.Nodegraph(6, 1e6, 2)
    try:
        readsparser = khmer.ReadParser(utils.get_test_data("test-empty.fa"))
        nodegraph.consume_seqfile_and_tag(readsparser)
        assert 0, "this should fail"
    except OSError as e:
        print(str(e))
    except ValueError as e:
        print(str(e))


def test_n_occupied_save_load():
    filename = utils.get_test_data('random-20-a.fa')

    nodegraph = khmer.Nodegraph(20, 100000, 3)

    for _, record in enumerate(screed.open(filename)):
        nodegraph.consume(record.sequence)

    assert nodegraph.n_occupied() == 3884
    assert nodegraph.n_unique_kmers() == 3960

    savefile = utils.get_temp_filename('out')
    nodegraph.save(savefile)

    ng2 = Nodegraph.load(savefile)
    assert ng2.n_occupied() == 3884, ng2.n_occupied()
    assert ng2.n_unique_kmers() == 0    # this is intended behavior, sigh.


def test_n_occupied_vs_countgraph():
    filename = utils.get_test_data('random-20-a.fa')

    nodegraph = khmer.Nodegraph(20, 100000, 3)
    countgraph = khmer.Countgraph(20, 100000, 3)

    assert nodegraph.n_occupied() == 0, nodegraph.n_occupied()
    assert countgraph.n_occupied() == 0, countgraph.n_occupied()

    assert nodegraph.n_unique_kmers() == 0, nodegraph.n_unique_kmers()
    assert countgraph.n_unique_kmers() == 0, countgraph.n_unique_kmers()

    for _, record in enumerate(screed.open(filename)):
        nodegraph.consume(record.sequence)
        countgraph.consume(record.sequence)

    assert nodegraph.hashsizes() == nodegraph.hashsizes()

    # these are all the same -- good :).
    assert nodegraph.n_occupied() == 3884, nodegraph.n_occupied()
    assert countgraph.n_occupied() == 3884, countgraph.n_occupied()

    assert nodegraph.n_unique_kmers() == 3960, nodegraph.n_unique_kmers()
    assert countgraph.n_unique_kmers() == 3960, countgraph.n_unique_kmers()


def test_n_occupied_vs_countgraph_another_size():
    filename = utils.get_test_data('random-20-a.fa')

    nodegraph = khmer.Nodegraph(20, 10000, 3)
    countgraph = khmer.Countgraph(20, 10000, 3)

    assert nodegraph.n_occupied() == 0, nodegraph.n_occupied()
    assert countgraph.n_occupied() == 0, countgraph.n_occupied()

    assert nodegraph.n_unique_kmers() == 0, nodegraph.n_unique_kmers()
    assert countgraph.n_unique_kmers() == 0, countgraph.n_unique_kmers()

    for _, record in enumerate(screed.open(filename)):
        nodegraph.consume(record.sequence)
        countgraph.consume(record.sequence)

    assert nodegraph.hashsizes() == nodegraph.hashsizes()

    # these are all the same -- good :).
    assert nodegraph.n_occupied() == 3269, nodegraph.n_occupied()
    assert countgraph.n_occupied() == 3269, countgraph.n_occupied()

    assert nodegraph.n_unique_kmers() == 3916, nodegraph.n_unique_kmers()
    assert countgraph.n_unique_kmers() == 3916, countgraph.n_unique_kmers()


def test_traverse_linear_path():
    contigfile = utils.get_test_data('simple-genome.fa')
    contig = list(screed.open(contigfile))[0].sequence

    K = 21

    nodegraph = khmer.Nodegraph(K, 1e5, 4)
    stopgraph = khmer.Nodegraph(K, 1e5, 4)

    nodegraph.consume(contig)

    degree_nodes = khmer.HashSet(K)
    size, conns, visited = nodegraph.traverse_linear_path(contig[:K],
                                                          degree_nodes,
                                                          stopgraph)
    assert size == 980
    assert len(conns) == 0
    assert len(visited) == 980


def test_find_high_degree_nodes():
    contigfile = utils.get_test_data('simple-genome.fa')
    contig = list(screed.open(contigfile))[0].sequence

    K = 21

    nodegraph = khmer.Nodegraph(K, 1e5, 4)
    stopgraph = khmer.Nodegraph(K, 1e5, 4)

    nodegraph.consume(contig)

    degree_nodes = nodegraph.find_high_degree_nodes(contig)
    assert len(degree_nodes) == 0


def test_find_high_degree_nodes_2():
    contigfile = utils.get_test_data('simple-genome.fa')
    contig = list(screed.open(contigfile))[0].sequence

    K = 21

    nodegraph = khmer.Nodegraph(K, 1e5, 4)

    nodegraph.consume(contig)
    nodegraph.count(contig[2:22] + 'G')   # will add another neighbor to 1:22
    print(nodegraph.neighbors(contig[1:22]))

    degree_nodes = nodegraph.find_high_degree_nodes(contig)
    assert len(degree_nodes) == 1
    assert nodegraph.hash(contig[1:22]) in degree_nodes


def test_traverse_linear_path_2():
    contigfile = utils.get_test_data('simple-genome.fa')
    contig = list(screed.open(contigfile))[0].sequence
    print('contig len', len(contig))

    K = 21

    nodegraph = khmer.Nodegraph(K, 1e5, 4)
    stopgraph = khmer.Nodegraph(K, 1e5, 4)

    nodegraph.consume(contig)
    nodegraph.count(contig[101:121] + 'G')  # will add another neighbor
    print(nodegraph.neighbors(contig[101:122]))

    degree_nodes = nodegraph.find_high_degree_nodes(contig)

    assert len(degree_nodes) == 1
    assert nodegraph.hash(contig[100:121]) in degree_nodes

    # traverse from start, should end at node 100:121
    size, conns, visited = nodegraph.traverse_linear_path(contig[0:21],
                                                          degree_nodes,
                                                          stopgraph)

    print(size, list(conns), list(visited))
    assert size == 100
    assert len(visited) == 100
    assert nodegraph.hash(contig[100:121]) in conns
    assert len(conns) == 1

    for k in conns:                       # everything in connections => stop
        assert stopgraph.get(k)

    for k in visited:                     # nothing in visited => stop
        assert not stopgraph.get(k)

    # traverse from immediately after 100:121, should end at the end
    size, conns, visited = nodegraph.traverse_linear_path(contig[101:122],
                                                          degree_nodes,
                                                          stopgraph)

    print(size, list(conns), list(visited))
    assert size == 879
    assert len(visited) == 879
    assert nodegraph.hash(contig[100:121]) in conns
    assert len(conns) == 1

    for k in conns:                       # everything in connections => stop
        assert stopgraph.get(k)

    for k in visited:                     # nothing in visited => stop
        assert not stopgraph.get(k)

    # traverse from end, should end at 100:121
    size, conns, visited = nodegraph.traverse_linear_path(contig[-21:],
                                                          degree_nodes,
                                                          stopgraph)

    print(size, list(conns), len(visited))
    assert size == 879
    assert len(visited) == 879
    assert nodegraph.hash(contig[100:121]) in conns
    assert len(conns) == 1

    for k in conns:                       # everything in connections => stop
        assert stopgraph.get(k)

    for k in visited:                     # nothing in visited => stop
        assert not stopgraph.get(k)


def test_traverse_linear_path_3_stopgraph():
    contigfile = utils.get_test_data('simple-genome.fa')
    contig = list(screed.open(contigfile))[0].sequence
    print('contig len', len(contig))

    K = 21

    nodegraph = khmer.Nodegraph(K, 1e5, 4)
    stopgraph = khmer.Nodegraph(K, 1e5, 4)

    nodegraph.consume(contig)
    nodegraph.count(contig[101:121] + 'G')  # will add another neighbor
    print(nodegraph.neighbors(contig[101:122]))

    degree_nodes = nodegraph.find_high_degree_nodes(contig)

    assert len(degree_nodes) == 1
    assert nodegraph.hash(contig[100:121]) in degree_nodes

    stopgraph.count(contig[101:122])       # stop traversal - only adj to start

    size, conns, visited = nodegraph.traverse_linear_path(contig[101:122],
                                                          degree_nodes,
                                                          stopgraph)

    print(size, list(conns), len(visited))
    assert size == 0
    assert len(visited) == 0
    assert len(conns) == 0


@pytest.mark.parametrize('ntables,targetsize', [
    (4, 1e5),
    (6, 1e5),
    (8, 1e5),
    (5, 1e6),
    (7, 1e6),
    (9, 1e6),
])
def test_create_matching_nodegraph(ntables, targetsize):
    cg = khmer.Countgraph(31, targetsize, ntables)
    ng = create_matching_nodegraph(cg)
    assert cg.hashsizes() == ng.hashsizes()
