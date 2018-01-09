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
# pylint: disable=missing-docstring,invalid-name,no-member,no-self-use
# pylint: disable=protected-access

import khmer
from khmer._oxli.legacy_partitioning import SubsetPartition, PrePartitionInfo
import screed

import os
from . import khmer_tst_utils as utils


def teardown():
    utils.cleanup()


class Test_RandomData(object):

    def test_3_merge_013(self):
        ht = khmer.Nodegraph(20, 4 ** 4 + 1, 2)

        filename = utils.get_test_data('test-graph2.fa')

        (total_reads, _) = ht.consume_seqfile_and_tag(filename)
        assert total_reads == 3, total_reads

        divvy = ht.divide_tags_into_subsets(1)
        assert len(divvy) is 3
        (a, b, _) = divvy

        x = ht.do_subset_partition(a, a)
        ht.merge_subset(x)

        y = ht.do_subset_partition(b, 0)
        ht.merge_subset(y)

        outfile = utils.get_temp_filename('out')
        n_partitions = ht.output_partitions(filename, outfile)
        assert n_partitions == 1, n_partitions        # combined.

    def test_3_merge_023(self):
        ht = khmer.Nodegraph(20, 4 ** 4 + 1, 2)
        filename = utils.get_test_data('test-graph2.fa')

        (total_reads, _) = ht.consume_seqfile_and_tag(filename)
        assert total_reads == 3, total_reads

        divvy = ht.divide_tags_into_subsets(1)
        assert len(divvy) is 3
        (a, b, c) = divvy

        x = ht.do_subset_partition(b, c)
        ht.merge_subset(x)

        y = ht.do_subset_partition(a, b)
        ht.merge_subset(y)

        outfile = utils.get_temp_filename('out.part')
        n_partitions = ht.output_partitions(filename, outfile)
        assert n_partitions == 1, n_partitions        # combined.

    def test_5_merge_046(self):
        ht = khmer.Nodegraph(20, 4 ** 4 + 1, 2)
        filename = utils.get_test_data('test-graph5.fa')

        (total_reads, _) = ht.consume_seqfile_and_tag(filename)
        assert total_reads == 6, total_reads

        divvy = ht.divide_tags_into_subsets(1)
        divvy = list(divvy)

        x = ht.do_subset_partition(divvy[0], divvy[4])
        ht.merge_subset(x)

        y = ht.do_subset_partition(divvy[4], 0)
        ht.merge_subset(y)

        outfile = utils.get_temp_filename('out.part')
        n_partitions = ht.output_partitions(filename, outfile)
        assert n_partitions == 1, n_partitions        # combined.

    def test_random_20_a_succ(self):
        ht = khmer.Nodegraph(20, 4 ** 7 + 1, 2)
        filename = utils.get_test_data('random-20-a.fa')
        outfile = utils.get_temp_filename('out')

        total_reads, _ = ht.consume_seqfile_and_tag(filename)

        subset_size = total_reads // 2 + total_reads % 2
        divvy = ht.divide_tags_into_subsets(subset_size)
        divvy = list(divvy)
        assert len(divvy) == 4

        x = ht.do_subset_partition(divvy[0], divvy[2])
        ht.merge_subset(x)
        y = ht.do_subset_partition(divvy[2], 0)
        ht.merge_subset(y)

        n_partitions = ht.output_partitions(filename, outfile)
        assert n_partitions == 1, n_partitions

    def test_random_20_a_succ_II(self):
        ht = khmer.Nodegraph(20, 4 ** 7 + 1, 2)
        filename = utils.get_test_data('random-20-a.fa')
        outfile = utils.get_temp_filename('out')

        total_reads, _ = ht.consume_seqfile_and_tag(filename)

        subset_size = total_reads // 2 + total_reads % 2
        divvy = ht.divide_tags_into_subsets(subset_size)
        divvy = list(divvy)
        assert len(divvy) == 4

        x = ht.do_subset_partition(divvy[0], divvy[2])
        y = ht.do_subset_partition(divvy[2], 0)
        ht.merge_subset(x)
        ht.merge_subset(y)

        n_partitions = ht.output_partitions(filename, outfile)
        assert n_partitions == 1, n_partitions

    def test_random_20_a_succ_III(self):
        ht = khmer.Nodegraph(20, 4 ** 7 + 1, 2)
        filename = utils.get_test_data('random-20-a.fa')
        outfile = utils.get_temp_filename('out')

        total_reads, _ = ht.consume_seqfile_and_tag(filename)

        subset_size = total_reads // 2 + total_reads % 2
        divvy = ht.divide_tags_into_subsets(subset_size)
        divvy = list(divvy)
        assert len(divvy) == 4, len(divvy)

        x = ht.do_subset_partition(divvy[0], divvy[2])
        y = ht.do_subset_partition(divvy[2], 0)

        x._validate_partitionmap()
        y._validate_partitionmap()

        ht.merge_subset(y)
        ht.merge_subset(x)

        n_partitions = ht.output_partitions(filename, outfile)
        assert n_partitions == 1, n_partitions

    def test_random_20_a_succ_IV(self):
        ht = khmer.Nodegraph(20, 4 ** 7 + 1, 2)
        filename = utils.get_test_data('random-20-a.fa')
        outfile = utils.get_temp_filename('out')

        ht.consume_seqfile_and_tag(filename)
        subsets = []

        divvy = ht.divide_tags_into_subsets(1)
        divvy = list(divvy)
        divvy.append(0)
        for i in range(len(divvy) - 1):
            x = ht.do_subset_partition(divvy[i], divvy[i + 1])
            subsets.append(x)

        for x in reversed(subsets):
            ht.merge_subset(x)

        n_partitions = ht.output_partitions(filename, outfile)
        assert n_partitions == 1, n_partitions

    def test_random_20_a_succ_IV_save(self):
        ht = khmer.Nodegraph(20, 4 ** 7 + 1, 2)
        filename = utils.get_test_data('random-20-a.fa')

        savefile_ht = utils.get_temp_filename('ht')
        savefile_tags = utils.get_temp_filename('tags')
        outfile = filename + utils.get_temp_filename('out')

        ht.consume_seqfile_and_tag(filename)

        ht.save(savefile_ht)
        ht.save_tagset(savefile_tags)

        del ht
        ht = khmer.Nodegraph(20, 4 ** 7 + 1, 2)

        ht = khmer.Nodegraph.load(savefile_ht)
        ht.load_tagset(savefile_tags)

        divvy = ht.divide_tags_into_subsets(1)
        divvy = list(divvy)
        divvy.append(0)

        subsets = []
        for i in range(len(divvy) - 1):
            x = ht.do_subset_partition(divvy[i], divvy[i + 1])
            subsets.append(x)

        for x in reversed(subsets):
            ht.merge_subset(x)

        n_partitions = ht.output_partitions(filename, outfile)
        assert n_partitions == 1, n_partitions


class Test_SaveLoadPmap(object):

    def test_save_load_merge(self):
        ht = khmer.Nodegraph(20, 4 ** 4 + 1, 2)
        filename = utils.get_test_data('test-graph2.fa')

        (total_reads, _) = ht.consume_seqfile_and_tag(filename)
        assert total_reads == 3, total_reads

        divvy = ht.divide_tags_into_subsets(1)
        assert len(divvy) == 3
        (a, b, _) = divvy

        outfile1 = utils.get_temp_filename('x.pmap')
        outfile2 = utils.get_temp_filename('y.pmap')

        x = ht.do_subset_partition(a, b)
        x.save_partitionmap(outfile1)
        del x

        y = ht.do_subset_partition(b, 0)
        y.save_partitionmap(outfile2)
        del y

        a = SubsetPartition.load(outfile1, ht)
        b = SubsetPartition.load(outfile2, ht)

        ht.merge_subset(a)
        ht.merge_subset(b)

        outfile = utils.get_temp_filename('out.part')
        n_partitions = ht.output_partitions(filename, outfile)
        assert n_partitions == 1, n_partitions        # combined.

    def test_save_load_merge_truncate(self):
        ht = khmer.Nodegraph(20, 4 ** 4 + 1, 2)
        filename = utils.get_test_data('test-graph2.fa')

        (total_reads, _) = ht.consume_seqfile_and_tag(filename)
        assert total_reads == 3, total_reads

        divvy = ht.divide_tags_into_subsets(1)
        print(divvy)
        assert len(divvy) is 3
        (a, b, _) = divvy

        outfile1 = utils.get_temp_filename('x.pmap')
        outfile2 = utils.get_temp_filename('y.pmap')

        x = ht.do_subset_partition(a, b)
        x.save_partitionmap(outfile1)
        del x

        y = ht.do_subset_partition(b, 0)
        y.save_partitionmap(outfile2)
        del y

        outfile3 = utils.get_temp_filename('z.pmap')
        data = open(outfile1, 'rb').read()

        for i in range(len(data)):
            fp = open(outfile3, 'wb')
            fp.write(data[:i])
            fp.close()

            try:
                a = SubsetPartition.load(outfile3, ht)
                assert 0, "this should not pass"
            except OSError as err:
                print(str(err), i)

    def test_save_load_merge_2(self):
        ht = khmer.Nodegraph(20, 4 ** 8 + 1, 2)
        filename = utils.get_test_data('random-20-a.fa')

        (total_reads, _) = ht.consume_seqfile_and_tag(filename)

        subset_size = total_reads // 2 + total_reads % 2
        divvy = ht.divide_tags_into_subsets(subset_size)
        divvy = list(divvy)

        outfile1 = utils.get_temp_filename('x.pmap')
        outfile2 = utils.get_temp_filename('y.pmap')

        x = ht.do_subset_partition(divvy[0], divvy[1])
        x.save_partitionmap(outfile1)
        del x

        y = ht.do_subset_partition(divvy[1], 0)
        y.save_partitionmap(outfile2)
        del y

        assert os.path.exists(outfile1)
        assert os.path.exists(outfile2)
        a = SubsetPartition.load(outfile1, ht)
        b = SubsetPartition.load(outfile2, ht)

        ht.merge_subset(a)
        ht.merge_subset(b)

        outfile = utils.get_temp_filename('out.part')
        n_partitions = ht.output_partitions(filename, outfile)
        assert n_partitions == 1, n_partitions        # combined.

    def test_save_load_merge_nexist(self):
        ht = khmer.Nodegraph(20, 1, 1)
        try:
            ht.load_partitionmap('this does not exist')
            assert 0, "this should not succeed"
        except OSError as e:
            print(str(e))

    def test_save_merge_from_disk(self):
        ht = khmer.Nodegraph(20, 4 ** 4 + 1, 2)
        filename = utils.get_test_data('test-graph2.fa')

        (total_reads, _) = ht.consume_seqfile_and_tag(filename)
        assert total_reads == 3, total_reads

        divvy = ht.divide_tags_into_subsets(1)
        print(divvy)
        (a, b, _) = divvy

        outfile1 = utils.get_temp_filename('x.pmap')
        outfile2 = utils.get_temp_filename('y.pmap')

        x = ht.do_subset_partition(a, b)
        x.save_partitionmap(outfile1)
        del x

        y = ht.do_subset_partition(b, 0)
        y.save_partitionmap(outfile2)
        del y

        ht.merge_subset_from_disk(outfile1)
        ht.merge_subset_from_disk(outfile2)

        outfile = utils.get_temp_filename('out.part')
        n_partitions = ht.output_partitions(filename, outfile)
        assert n_partitions == 1, n_partitions        # combined.

    def test_save_merge_from_disk_2(self):
        ht = khmer.Nodegraph(20, 4 ** 7 + 1, 2)
        filename = utils.get_test_data('random-20-a.fa')

        (total_reads, _) = ht.consume_seqfile_and_tag(filename)

        subset_size = total_reads // 2 + total_reads % 2
        divvy = ht.divide_tags_into_subsets(subset_size)
        divvy = list(divvy)

        outfile1 = utils.get_temp_filename('x.pmap')
        outfile2 = utils.get_temp_filename('y.pmap')

        x = ht.do_subset_partition(divvy[0], divvy[1])
        x.save_partitionmap(outfile1)
        del x

        y = ht.do_subset_partition(divvy[1], 0)
        y.save_partitionmap(outfile2)
        del y

        assert os.path.exists(outfile1)
        assert os.path.exists(outfile2)
        ht.merge_subset_from_disk(outfile1)
        ht.merge_subset_from_disk(outfile2)

        outfile = utils.get_temp_filename('out.part')
        n_partitions = ht.output_partitions(filename, outfile)
        assert n_partitions == 1, n_partitions        # combined.

    def test_save_merge_from_disk_file_not_exist(self):
        ht = khmer.Nodegraph(20, 4 ** 4 + 1, 2)
        filename = utils.get_test_data('test-graph2.fa')

        (total_reads, _) = ht.consume_seqfile_and_tag(filename)
        assert total_reads == 3, total_reads

        outfile1 = utils.get_temp_filename('x.pmap')

        # fail to create file... => failure expected

        try:
            ht.merge_subset_from_disk(outfile1)
            assert 0, "this should fail"
        except OSError as e:
            print(str(e))

    def test_merge_from_disk_file_bad_type(self):
        ht = khmer.Nodegraph(20, 4 ** 4 + 1, 2)
        infile = utils.get_test_data('goodversion-k12.ht')

        try:
            ht.merge_subset_from_disk(infile)
            assert 0, "this should fail"
        except OSError as e:
            print(str(e))

    def test_merge_from_disk_file_version(self):
        ht = khmer.Nodegraph(20, 4 ** 4 + 1, 2)
        infile = utils.get_test_data('badversion-k12.ht')

        try:
            ht.merge_subset_from_disk(infile)
            assert 0, "this should fail"
        except OSError as e:
            print(str(e))

    def test_save_merge_from_disk_ksize(self):
        ht = khmer.Nodegraph(20, 4 ** 4 + 1, 2)
        filename = utils.get_test_data('test-graph2.fa')

        (total_reads, _) = ht.consume_seqfile_and_tag(filename)
        assert total_reads == 3, total_reads

        divvy = ht.divide_tags_into_subsets(1)
        print(divvy)
        (a, b, _) = divvy

        outfile1 = utils.get_temp_filename('x.pmap')
        x = ht.do_subset_partition(a, b)
        x.save_partitionmap(outfile1)
        del x

        ht = khmer.Nodegraph(19, 1, 1)
        try:
            ht.merge_subset_from_disk(outfile1)
            assert 0, "this should fail"
        except OSError as e:
            print(str(e))


def test_save_load_merge_on_graph():
    ht = khmer.Nodegraph(20, 4 ** 4 + 1, 2)
    filename = utils.get_test_data('test-graph2.fa')

    (total_reads, _) = ht.consume_seqfile_and_tag(filename)
    assert total_reads == 3, total_reads

    divvy = ht.divide_tags_into_subsets(1)
    print(divvy)
    assert len(divvy) is 3
    (a, b, _) = divvy

    outfile1 = utils.get_temp_filename('x.pmap')
    outfile2 = utils.get_temp_filename('y.pmap')

    x = ht.do_subset_partition(a, b)
    x.save_partitionmap(outfile1)
    del x

    y = ht.do_subset_partition(b, 0)
    y.save_partitionmap(outfile2)
    del y

    a = ht.load_partitionmap(outfile1)  # <-- this is different
    b = SubsetPartition.load(outfile2, ht)

    ht.merge_subset(b)

    outfile = utils.get_temp_filename('out.part')
    n_partitions = ht.output_partitions(filename, outfile)
    assert n_partitions == 1, n_partitions        # combined.


def test_save_load_on_graph_truncate():
    ht = khmer.Nodegraph(20, 4 ** 4 + 1, 2)
    filename = utils.get_test_data('test-graph2.fa')

    (total_reads, _) = ht.consume_seqfile_and_tag(filename)
    assert total_reads == 3, total_reads

    divvy = ht.divide_tags_into_subsets(1)
    print(divvy)
    assert len(divvy) is 3
    (a, b, _) = divvy

    outfile1 = utils.get_temp_filename('x.pmap')
    outfile2 = utils.get_temp_filename('y.pmap')

    x = ht.do_subset_partition(a, b)
    x.save_partitionmap(outfile1)
    del x

    y = ht.do_subset_partition(b, 0)
    y.save_partitionmap(outfile2)
    del y

    outfile3 = utils.get_temp_filename('z.pmap')
    data = open(outfile1, 'rb').read()

    for i in range(len(data)):
        fp = open(outfile3, 'wb')
        fp.write(data[:i])
        fp.close()

        try:
            a = ht.load_partitionmap(outfile3)
            assert 0, "this should not pass"
        except OSError as err:
            print(str(err), i)


def test_output_partitions():
    filename = utils.get_test_data('test-output-partitions.fa')

    ht = khmer.Nodegraph(10, 1, 1)
    ht.set_partition_id('TTAGGACTGC', 2)
    ht.set_partition_id('TGCGTTTCAA', 3)
    ht.set_partition_id('ATACTGTAAA', 4)

    outfile = utils.get_temp_filename('part')
    ht.output_partitions(filename, outfile)

    data = open(outfile).read()
    assert len(data)

    records = [r for r in screed.open(outfile)]
    names = [r.name for r in records]
    parts = [n.rsplit('\t', 1)[1] for n in names]

    assert parts[0] == '2'
    assert parts[1] == '3'
    assert parts[2] == '4'


test_output_partitions.runme = True


def test_tiny_real_partitions():
    filename = utils.get_test_data('real-partition-tiny.fa')

    ht = khmer.Nodegraph(32, 8e2, 4)
    ht.consume_seqfile_and_tag(filename)

    subset = ht.do_subset_partition(0, 0)
    ht.merge_subset(subset)

    outfile = utils.get_temp_filename('part')
    ht.output_partitions(filename, outfile)

    data = open(outfile).read()

    assert len(data)

    records = [r for r in screed.open(outfile)]
    names = [r.name for r in records]
    parts = [n.rsplit('\t', 1)[1] for n in names]

    assert len(parts) == 2, len(parts)
    assert len(set(parts)) == 1
    assert set(parts) != set(['0'])

    test_tiny_real_partitions.runme = True


def test_small_real_partitions():
    filename = utils.get_test_data('real-partition-small.fa')

    ht = khmer.Nodegraph(32, 2e3, 4)
    ht.consume_seqfile_and_tag(filename)

    subset = ht.do_subset_partition(0, 0)
    ht.merge_subset(subset)

    outfile = utils.get_temp_filename('part')
    ht.output_partitions(filename, outfile)

    data = open(outfile).read()
    assert len(data)

    records = [r for r in screed.open(outfile)]
    names = [r.name for r in records]
    parts = [n.rsplit('\t', 1)[1] for n in names]

    assert len(parts) == 6, len(parts)
    assert len(set(parts)) == 1
    assert set(parts) != set(['0'])


test_small_real_partitions.runme = True

first = """\
CAGACTTGGAAGCTGAGAGTCCGACGTCACTGCCTCAACTCGCGCAAATGTTCCCGCCAA\
ATTGTATCCTAGGGATCTTCCATAAGCTTATATACGGGGGTTTCCAAGGCCCTGATGCCA\
GTGCCTAATCTTTTGGAGTCCTCTCAGGGCCACTAGATGCCATGCTACGCGTCCCAGGTT\
GGCCTGAGGGTCTACACGGAGTGGGAAGCATGGGTACCTTAGCGAACATTCATACTGGCC\
TGTTTATGCTTATCAGACTTCAGCTTCGCTTAGCGCGTCACCGTTTGTAACTTGTTATCT\
"""

second = """\
TGTTTATGCTTATCAGACTTCAGCTTCGCTTAGCGCGTCACCGTTTGTAACTTGTTATCT\
GACTGTAGACTTGAACCTCGATGGAATGCAGGTCCCATTCTCTGGCCTGACTCATGGAAC\
CGAGGCCAAAAAAGCATGGCACGAAGACGCTATGCGAGGGTGCTCGCCCATGTCGTCGCC\
GTACCACGACAGATTTATACAATGCGTTTCTACAGGCCCCATTGGGAACAAACAAAAAGT\
CCTCGGGCCTTTCCGTTCCGTTGCCGCCCAAGCTCTCTAGCATCGAATCGGTCAAGCGGT\
"""


def test_partition_on_abundance_1():
    kh = khmer.Countgraph(20, 1e3, 4)
    for _ in range(10):
        print(kh.consume_and_tag(first))

    for _ in range(10):
        print(kh.consume_and_tag(second))

    # all paths in 'a' and 'b'
    p = kh.do_subset_partition_with_abundance(10, 50)
    x = p.count_partitions()
    assert x == (1, 0)                  # one partition, no remainders


def test_partition_on_abundance_2():
    kh = khmer.Countgraph(20, 1e3, 4)
    for _ in range(10):
        print(kh.consume_and_tag(first))

    for _ in range(5):
        print(kh.consume_and_tag(second))

    # all paths in 'a'
    p = kh.do_subset_partition_with_abundance(10, 50)
    x = p.count_partitions()
    assert x == (1, 6)                  # one partition, six disconnected


def test_partition_on_abundance_3():
    kh = khmer.Countgraph(20, 1e4, 4)
    for _ in range(10):
        print(kh.consume_and_tag(first))

    for _ in range(5):
        print(kh.consume_and_tag(second))

    # this will get paths only in 'a'
    p = kh.do_subset_partition_with_abundance(10, 50)

    # this will get paths only in 'b'
    p = kh.do_subset_partition_with_abundance(5, 10)

    x = p.count_partitions()
    print(x)
    assert x == (2, 2)                  # two partitions, two ignored tags


def test_partition_overlap_2():
    kh = khmer.Countgraph(20, 1e4, 4)
    for _ in range(10):
        kh.consume_and_tag(first)

    for _ in range(5):
        kh.consume_and_tag(second)

    # this will get paths only in 'a'
    p1 = kh.do_subset_partition_with_abundance(10, 50)

    # this will get paths only in 'b'
    p2 = kh.do_subset_partition_with_abundance(5, 10)

    # p1.report_on_partitions()
    # p2.report_on_partitions()

    x = p1.partition_sizes()
    assert x == ([(3, 8)], 0), x

    x = p2.partition_sizes()
    x[0].sort(key=lambda pair: pair[0])
    assert x == ([(3, 6), (5, 6)], 0), x

    x = p1.partition_average_coverages(kh)
    assert x == [(3, 11)]

    x = p2.partition_average_coverages(kh)
    x.sort(key=lambda pair: pair[0])
    assert x == [(3, 5), (5, 10)], x
