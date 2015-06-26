from __future__ import print_function
from __future__ import absolute_import
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
import khmer
import screed

from . import khmer_tst_utils as utils

from nose.plugins.attrib import attr


def teardown():
    utils.cleanup()


class Test_ExactGraphFu(object):

    def setup(self):
        self.ht = khmer.new_hashbits(12, 1e4)

    def test_counts(self):
        ht = self.ht
        ht.consume_fasta(utils.get_test_data('test-graph.fa'))

        kmer = "TTAGGACTGCAC"
        x = ht.calc_connected_graph_size(kmer)
        assert x == 69, x

        kmer = "TGCGTTTCAATC"
        x = ht.calc_connected_graph_size(kmer)
        assert x == 68, x

        kmer = "ATACTGTAAATA"
        x = ht.calc_connected_graph_size(kmer)
        assert x == 36, x

    def test_graph_links_next_a(self):
        ht = self.ht
        word = "TGCGTTTCAATC"
        ht.consume(word)
        ht.consume(word[1:] + "A")

        x = ht.calc_connected_graph_size(word)
        assert x == 2

    def test_graph_links_next_c(self):
        ht = self.ht
        word = "TGCGTTTCAATC"
        ht.consume(word)
        ht.consume(word[1:] + "C")

        x = ht.calc_connected_graph_size(word)
        assert x == 2

    def test_graph_links_next_g(self):
        ht = self.ht
        word = "TGCGTTTCAATC"
        ht.consume(word)
        ht.consume(word[1:] + "G")

        x = ht.calc_connected_graph_size(word)
        assert x == 2

    def test_graph_links_next_t(self):
        ht = self.ht
        word = "TGCGTTTCAATC"
        ht.consume(word)
        ht.consume(word[1:] + "T")

        x = ht.calc_connected_graph_size(word)
        assert x == 2

    def test_graph_links_prev_a(self):
        ht = self.ht
        word = "TGCGTTTCAATC"
        ht.consume(word)
        ht.consume("A" + word[:-1])

        x = ht.calc_connected_graph_size(word)
        assert x == 2

    def test_graph_links_prev_c(self):
        ht = self.ht
        word = "TGCGTTTCAATC"
        ht.consume(word)
        ht.consume("C" + word[:-1])

        x = ht.calc_connected_graph_size(word)
        assert x == 2

    def test_graph_links_prev_g(self):
        ht = self.ht
        word = "TGCGTTTCAATC"
        ht.consume(word)
        ht.consume("G" + word[:-1])

        x = ht.calc_connected_graph_size(word)
        assert x == 2

    def test_graph_links_prev_t(self):
        ht = self.ht
        word = "TGCGTTTCAATC"
        ht.consume(word)
        ht.consume("T" + word[:-1])

        x = ht.calc_connected_graph_size(word)
        assert x == 2


class Test_InexactGraphFu(object):

    def setup(self):
        self.ht = khmer.new_hashbits(12, 4 ** 2 + 1)

    def test_graph_links_next_a(self):
        ht = self.ht
        word = "TGCGTTTCAATC"
        ht.consume(word)
        ht.consume(word[1:] + "A")

        x = ht.calc_connected_graph_size(word)
        assert x == 2

    def test_graph_links_next_c(self):
        ht = self.ht
        word = "TGCGTTTCAATC"
        ht.consume(word)
        ht.consume(word[1:] + "C")

        x = ht.calc_connected_graph_size(word)
        assert x == 2

    def test_graph_links_next_g(self):
        ht = self.ht
        word = "TGCGTTTCAATC"
        ht.consume(word)
        ht.consume(word[1:] + "G")

        x = ht.calc_connected_graph_size(word)
        assert x == 2

    def test_graph_links_next_t(self):
        ht = self.ht
        word = "TGCGTTTCAATC"
        ht.consume(word)
        ht.consume(word[1:] + "T")

        x = ht.calc_connected_graph_size(word)
        assert x == 2

    def test_graph_links_prev_a(self):
        ht = self.ht
        word = "TGCGTTTCAATC"
        ht.consume(word)
        ht.consume("A" + word[:-1])

        x = ht.calc_connected_graph_size(word)
        assert x == 2

    def test_graph_links_prev_c(self):
        ht = self.ht
        word = "TGCGTTTCAATC"
        ht.consume(word)
        ht.consume("C" + word[:-1])

        x = ht.calc_connected_graph_size(word)
        assert x == 2

    def test_graph_links_prev_g(self):
        ht = self.ht
        word = "TGCGTTTCAATC"
        ht.consume(word)
        ht.consume("G" + word[:-1])

        x = ht.calc_connected_graph_size(word)
        assert x == 2

    def test_graph_links_prev_t(self):
        ht = self.ht
        word = "TGCGTTTCAATC"
        ht.consume(word)
        ht.consume("T" + word[:-1])

        x = ht.calc_connected_graph_size(word)
        assert x == 2

#


class Test_Partitioning(object):

    def test_output_unassigned(self):
        import screed

        filename = utils.get_test_data('random-20-a.fa')

        ht = khmer.new_hashbits(21, 4, 4)
        ht.consume_fasta_and_tag(filename)

        output_file = utils.get_temp_filename('part0test')
        ht.output_partitions(filename, output_file, True)

        len1 = len(list(screed.open(filename)))
        len2 = len(list(screed.open(output_file)))

        assert len1 > 0
        assert len1 == len2, (len1, len2)

    def test_not_output_unassigned(self):
        import screed

        filename = utils.get_test_data('random-20-a.fa')

        ht = khmer.new_hashbits(21, 4, 4)
        ht.consume_fasta_and_tag(filename)

        output_file = utils.get_temp_filename('parttest')
        ht.output_partitions(filename, output_file, False)

        len1 = len(list(screed.open(filename)))
        len2 = len(list(screed.open(output_file)))

        assert len1 > 0
        assert len2 == 0, len2

    def test_output_fq(self):
        filename = utils.get_test_data('random-20-a.fq')

        ht = khmer.new_hashbits(20, 1e4, 4)
        ht.consume_fasta_and_tag(filename)
        subset = ht.do_subset_partition(0, 0)
        ht.merge_subset(subset)

        output_file = utils.get_temp_filename('parttest')
        ht.output_partitions(filename, output_file, False)

        print(open(output_file).read())

        x = set([r.quality for r in screed.open(output_file)])
        assert x, x

    def test_disconnected_20_a(self):
        filename = utils.get_test_data('random-20-a.fa')

        ht = khmer.new_hashbits(21, 1e5, 4)
        ht.consume_fasta_and_tag(filename)

        subset = ht.do_subset_partition(0, 0)
        x = ht.subset_count_partitions(subset)
        assert x == (99, 0), x             # disconnected @ 21

    def test_connected_20_a(self):
        filename = utils.get_test_data('random-20-a.fa')

        ht = khmer.new_hashbits(20, 1e4, 4)
        ht.consume_fasta_and_tag(filename)

        subset = ht.do_subset_partition(0, 0)
        x = ht.subset_count_partitions(subset)
        assert x == (1, 0)             # connected @ 20

    def test_disconnected_20_b(self):
        filename = utils.get_test_data('random-20-b.fa')

        ht = khmer.new_hashbits(21, 1e4, 4)
        ht.consume_fasta_and_tag(filename)

        subset = ht.do_subset_partition(0, 0)
        x = ht.subset_count_partitions(subset)
        assert x == (99, 0), x             # disconnected @ 21

    def test_connected_20_b(self):
        filename = utils.get_test_data('random-20-b.fa')

        ht = khmer.new_hashbits(20, 1e4, 4)
        ht.consume_fasta_and_tag(filename)

        subset = ht.do_subset_partition(0, 0)
        x = ht.subset_count_partitions(subset)
        assert x == (1, 0)             # connected @ 20

    def test_disconnected_31_c(self):
        filename = utils.get_test_data('random-31-c.fa')

        ht = khmer.new_hashbits(32, 1e6, 4)
        ht.consume_fasta_and_tag(filename)

        subset = ht.do_subset_partition(0, 0)
        x = ht.subset_count_partitions(subset)
        assert x == (999, 0), x            # disconnected @ K = 32

    def test_connected_31_c(self):
        filename = utils.get_test_data('random-31-c.fa')

        ht = khmer.new_hashbits(31, 1e5, 4)
        ht.consume_fasta_and_tag(filename)

        subset = ht.do_subset_partition(0, 0)
        x = ht.subset_count_partitions(subset)
        assert x == (1, 0)             # connected @ K = 31

#


class Test_PythonAPI(object):

    def test_find_all_tags_kmersize(self):
        ht = khmer.new_hashbits(20, 4 ** 4 + 1)

        a = "ATTGGGACTCTGGGAGCACTTATCATGGAGAT"
        b = "GAGCACTTTAACCCTGCAGAGTGGCCAAGGCT"
        c = "GGAGCACTTATCATGGAGATATATCCCGTGCTTAAACATCGCACTTTAACCCTGCAGAGT"

        print(ht.consume(a))
        try:
            ppi = ht.find_all_tags(c[:19])
            assert False, "should raise a ValueError for wrong k-mer size"
        except ValueError:
            pass

        try:
            ppi = ht.find_all_tags(c[:21])
            assert False, "should raise a ValueError for wrong k-mer size"
        except ValueError:
            pass

    def test_ordered_connect(self):
        ht = khmer.new_hashbits(20, 4 ** 4 + 1)

        a = "ATTGGGACTCTGGGAGCACTTATCATGGAGAT"
        b = "GAGCACTTTAACCCTGCAGAGTGGCCAAGGCT"
        c = "GGAGCACTTATCATGGAGATATATCCCGTGCTTAAACATCGCACTTTAACCCTGCAGAGT"

        print(ht.consume(a))
        ppi = ht.find_all_tags(a[:20])
        pid = ht.assign_partition_id(ppi)
        assert pid == 0, pid

        print(ht.consume(b))
        ppi = ht.find_all_tags(b[:20])
        pid = ht.assign_partition_id(ppi)
        assert pid == 0, pid

        print(ht.consume(c))
        ppi = ht.find_all_tags(c[:20])
        pid = ht.assign_partition_id(ppi)
        assert pid == 2, pid

#
