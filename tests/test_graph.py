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

    def test_trim(self):
        ht = self.ht
        
        filename = os.path.join(thisdir, 'test-graph.fa')
        outfile = os.path.join(thisdir, 'test-graph.fa.out') # @CTB use tempfile
        ht.consume_fasta(filename)
        ht.trim_graphs(filename, 40, outfile)

        ht = khmer.new_hashbits(12, 4**12)
        ht.consume_fasta(outfile)

        x = ht.calc_connected_graph_size("TTAGGACTGCAC")
        assert x == 69, x
        
        x = ht.calc_connected_graph_size("TGCGTTTCAATC")
        assert x == 68, x
        
        x = ht.calc_connected_graph_size("ATACTGTAAATA")
        assert x == 0, x

    def test_graphsize_distrib(self):
        ht = self.ht
        ht.consume_fasta(os.path.join(thisdir, 'test-graph.fa'))
        x = ht.graphsize_distribution(200)

        assert sum(x) == 3, x
        assert x[69] == 1
        assert x[68] == 1
        assert x[36] == 1

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

    def test_trim(self):
        ht = self.ht
        filename = os.path.join(thisdir, 'test-graph.fa')
        outfile = os.path.join(thisdir, 'test-graph.fa.out')
        
        ht.consume_fasta(filename)
        ht.trim_graphs(filename, 40, outfile)

        ht = khmer.new_hashbits(12, 4**12)
        ht.consume_fasta(outfile)

        x = ht.calc_connected_graph_size("TTAGGACTGCAC")
        assert x >= 69, x
        
        x = ht.calc_connected_graph_size("TGCGTTTCAATC")
        assert x >= 68, x               # @CTB why 69??
        
        x = ht.calc_connected_graph_size("ATACTGTAAATA")
        assert x == 0, x

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

