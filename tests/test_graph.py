import os
thisdir = os.path.dirname(__file__)
thisdir = os.path.abspath(thisdir)

import khmer

class Test_ExactGraphFu(object):
    def setup(self):
        self.ht = khmer.new_hashbits(12, 4**12)

    def test_counts(self):
        ht = self.ht
        ht.consume_fasta(os.path.join(thisdir, 'test-graph.fa'))

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
        self.ht = khmer.new_hashbits(12, 4**8+1)

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

###

class Test_Partitioning(object):
    def test_disconnected_20_a(self):
        filename = os.path.join(thisdir, 'test-data/random-20-a.fa')

        ht = khmer.new_hashbits(21, 1e6, 4)
        ht.consume_fasta_and_tag(filename)

        subset = ht.do_subset_partition(0, 0)
        x = ht.subset_count_partitions(subset)
        assert x == (99, 0)             # disconnected @ 21

    def test_connected_20_a(self):
        filename = os.path.join(thisdir, 'test-data/random-20-a.fa')

        ht = khmer.new_hashbits(20, 1e6, 4)
        ht.consume_fasta_and_tag(filename)

        subset = ht.do_subset_partition(0, 0)
        x = ht.subset_count_partitions(subset)
        assert x == (1, 0)             # connected @ 20

    def test_disconnected_20_b(self):
        filename = os.path.join(thisdir, 'test-data/random-20-b.fa')

        ht = khmer.new_hashbits(21, 1e6, 4)
        ht.consume_fasta_and_tag(filename)

        subset = ht.do_subset_partition(0, 0)
        x = ht.subset_count_partitions(subset)
        assert x == (99, 0)             # disconnected @ 21

    def test_connected_20_b(self):
        filename = os.path.join(thisdir, 'test-data/random-20-b.fa')

        ht = khmer.new_hashbits(20, 1e6, 4)
        ht.consume_fasta_and_tag(filename)

        subset = ht.do_subset_partition(0, 0)
        x = ht.subset_count_partitions(subset)
        assert x == (1, 0)             # connected @ 20

    def test_disconnected_31_c(self):
        filename = os.path.join(thisdir, 'test-data/random-31-c.fa')

        ht = khmer.new_hashbits(32, 1e6, 4)
        ht.consume_fasta_and_tag(filename)

        subset = ht.do_subset_partition(0, 0)
        x = ht.subset_count_partitions(subset)
        assert x == (999, 0)            # disconnected @ K = 32

    def test_connected_31_c(self):
        filename = os.path.join(thisdir, 'test-data/random-31-c.fa')

        ht = khmer.new_hashbits(31, 1e6, 4)
        ht.consume_fasta_and_tag(filename)

        subset = ht.do_subset_partition(0, 0)
        x = ht.subset_count_partitions(subset)
        assert x == (1, 0)             # connected @ K = 31

###

class Test_PythonAPI(object):
    def test_ordered_connect(self):
        ht = khmer.new_hashbits(20, 4**15+1)

        a = "ATTGGGACTCTGGGAGCACTTATCATGGAGAT"
        b = "GAGCACTTTAACCCTGCAGAGTGGCCAAGGCT"
        c = "GGAGCACTTATCATGGAGATATATCCCGTGCTTAAACATCGCACTTTAACCCTGCAGAGT"

        print ht.consume(a)
        ppi = ht.find_all_tags(a[:20])
        pid = ht.assign_partition_id(ppi)
        assert pid == 0, pid
        
        print ht.consume(b)
        ppi = ht.find_all_tags(b[:20])
        pid = ht.assign_partition_id(ppi)
        assert pid == 0, pid
        
        print ht.consume(c)
        ppi = ht.find_all_tags(c[:20])
        pid = ht.assign_partition_id(ppi)
        assert pid == 2, pid

###

