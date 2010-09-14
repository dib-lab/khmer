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

class Test_PartitionCount(object):
    def test_simple_20_12(self):
        ht = khmer.new_hashbits(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph2.fa')
        ht.do_truncated_partition(filename, filename + '.out')
        n_partitions, n_unassigned, n_surrendered = ht.count_partitions()

        assert n_partitions == 1
        assert n_unassigned == 0
        assert n_surrendered == 0
        
    def test_simple_30_12(self):
        ht = khmer.new_hashbits(32, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph2.fa')
        ht.do_truncated_partition(filename, filename + '.out')
        n_partitions, n_unassigned, n_surrendered = ht.count_partitions()

        assert n_partitions == 0
        assert n_unassigned == 3
        assert n_surrendered == 0

    def test_surrendered(self):
        return                          # @@CTB
        ht = khmer.new_hashbits(32, 4**15+1)

        filename = os.path.join(thisdir, '../data/100k-surrendered.fa')
        ht.do_truncated_partition(filename, filename + '.out')
        n_partitions, n_unassigned, n_surrendered = ht.count_partitions()

        assert n_partitions == 1, n_partitions
        assert n_unassigned == 2, n_unassigned
        assert n_surrendered == 15, n_surrendered

### do_truncated_partition

class Test_SimpleConnectMe4(object):
    def test_simple_20_12(self):
        ht = khmer.new_hashbits(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph2.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out') # @CTB use tempfile
        n = ht.do_truncated_partition(filename, outfile)
        assert n == 1, n
        
    def test_simple_30_12(self):
        ht = khmer.new_hashbits(32, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph2.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile)
        assert n == 3, n

class Test_NoConnectMe4(object):
    def test_merge_20_12(self):
        ht = khmer.new_hashbits(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph3.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile)
        assert n == 2, n
        
    def test_merge_32_12(self):
        ht = khmer.new_hashbits(32, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph3.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile)
        assert n == 4, n

class Test_AnotherConnectMe4(object):
    def test_complex_20_12(self):
        ht = khmer.new_hashbits(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph4.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile)
        assert n == 2, n

    def test_complex_31_12(self):
        ht = khmer.new_hashbits(31, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph4.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile)
        assert n == 4, n

    def test_complex_32_12(self):
        ht = khmer.new_hashbits(32, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph4.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile)
        assert n == 5, n


class Test_MoreConnectMe4(object):
    def test_complex5_20_12(self):
        ht = khmer.new_hashbits(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph5.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile)
        assert n == 1, n

    def test_complex5_24_12(self):
        ht = khmer.new_hashbits(30, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph5.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile)
        assert n == 6, n

    def test_complex6_32_12(self):
        ht = khmer.new_hashbits(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph6.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile)
        assert n == 103, n

    def test_complex6_32_12_save(self):
        # this succeeds if save/load are null, or implemented. oops. @@CTB.
        ht = khmer.new_hashbits(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph6.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile)
        assert n == 103, n
        ht._validate_partitionmap()

        o1 = os.path.join(thisdir, 'xx.pmap')
        ht.save_partitionmap(o1)
        ht.load_partitionmap(o1)

        n = ht.do_truncated_partition(filename, outfile)
        assert n == 103, n
        ht._validate_partitionmap()

class Test_ThreadedSimpleConnectMe4(object):
    def test_simple_20_12(self):
        ht = khmer.new_hashbits(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph2.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out') # @CTB use tempfile
        n = ht.do_threaded_partition(filename)
        subset = ht.do_subset_partition(0, 0)
        n = ht.subset_count_partitions(subset)[0]
        assert n == 1, n
        
    def test_simple_30_12(self):
        ht = khmer.new_hashbits(32, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph2.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_threaded_partition(filename)
        subset = ht.do_subset_partition(0, 0)
        n = ht.subset_count_partitions(subset)[0]
        assert n == 0, n

class Test_ThreadedNoConnectMe4(object):
    def test_merge_20_12(self):
        ht = khmer.new_hashbits(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph3.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        ht.do_threaded_partition(filename)
        subset = ht.do_subset_partition(0, 0)
        n = ht.subset_count_partitions(subset)
        assert n == (1, 1, 0), n
        
    def test_merge_32_12(self):
        ht = khmer.new_hashbits(32, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph3.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        ht.do_threaded_partition(filename)
        subset = ht.do_subset_partition(0, 0)
        n = ht.subset_count_partitions(subset)
        assert n == (0, 4, 0), n

class Test_ThreadedAnotherConnectMe4(object):
    def test_complex_20_12(self):
        ht = khmer.new_hashbits(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph4.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        ht.do_threaded_partition(filename)
        subset = ht.do_subset_partition(0, 0)
        n = ht.subset_count_partitions(subset)[0]
        assert n == 2, n

    def test_complex_31_12(self):
        ht = khmer.new_hashbits(31, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph4.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        ht.do_threaded_partition(filename)
        subset = ht.do_subset_partition(0, 0)
        n = sum(ht.subset_count_partitions(subset))
        assert n == 4, n

    def test_complex_32_12(self):
        ht = khmer.new_hashbits(32, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph4.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        ht.do_threaded_partition(filename)
        subset = ht.do_subset_partition(0, 0)
        n = sum(ht.subset_count_partitions(subset))
        assert n == 5, n


class Test_ThreadedMoreConnectMe4(object):
    def test_complex5_20_12(self):
        ht = khmer.new_hashbits(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph5.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        ht.do_threaded_partition(filename)
        subset = ht.do_subset_partition(0, 0)
        n = ht.subset_count_partitions(subset)
        assert n == (1, 0, 0), n

    def test_complex5_24_12(self):
        ht = khmer.new_hashbits(30, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph5.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        ht.do_threaded_partition(filename)
        subset = ht.do_subset_partition(0, 0)
        n = ht.subset_count_partitions(subset)
        assert n == (0, 6, 0), n

    def test_complex6_32_12(self):
        ht = khmer.new_hashbits(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph6.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        ht.do_threaded_partition(filename)
        subset = ht.do_subset_partition(0, 0)
        n = ht.subset_count_partitions(subset)
        assert n == (102, 1, 0), n      # @@CTB?

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

class Test_RandomData(object):
    def test_random_20_a_succ(self):
        ht = khmer.new_hashbits(20, 4**13+1)
        filename = os.path.join(thisdir, 'test-data/random-20-a.fa')
        outfile = filename + '.out'

        n = ht.do_truncated_partition(filename, outfile)
        assert n == 1, n

    def test_random_20_a_fail(self):
        ht = khmer.new_hashbits(21, 4**13+1)
        filename = os.path.join(thisdir, 'test-data/random-20-a.fa')
        outfile = filename + '.out'

        n = ht.do_truncated_partition(filename, outfile)
        assert n == 99, n

    def test_random_20_b_succ(self):
        ht = khmer.new_hashbits(20, 4**13+1)
        filename = os.path.join(thisdir, 'test-data/random-20-b.fa')
        outfile = filename + '.out'

        n = ht.do_truncated_partition(filename, outfile)
        assert n == 1, n

    def test_random_20_b_fail(self):
        ht = khmer.new_hashbits(21, 4**13+1)
        filename = os.path.join(thisdir, 'test-data/random-20-b.fa')
        outfile = filename + '.out'

        n = ht.do_truncated_partition(filename, outfile)
        assert n == 99, n

    def test_random_31_a_succ(self):
        ht = khmer.new_hashbits(31, 4**14+1)
        filename = os.path.join(thisdir, 'test-data/random-31-c.fa')
        outfile = filename + '.out'

        n = ht.do_truncated_partition(filename, outfile)
        assert n == 1, n

    def test_random_31_b_fail(self):
        ht = khmer.new_hashbits(32, 4**14+1)
        filename = os.path.join(thisdir, 'test-data/random-31-c.fa')
        outfile = filename + '.out'

        n = ht.do_truncated_partition(filename, outfile)
        assert n == 999, n

class Test_Threaded_RandomData(object):
    def test_random_20_a_succ(self):
        ht = khmer.new_hashbits(20, 4**13+1)
        filename = os.path.join(thisdir, 'test-data/random-20-a.fa')
        outfile = filename + '.out'

        ht.do_threaded_partition(filename)
        subset = ht.do_subset_partition(0, 0)
        n = ht.subset_count_partitions(subset)
        assert n == (1, 0, 0), n

    def test_random_20_a_fail(self):
        ht = khmer.new_hashbits(21, 4**13+1)
        filename = os.path.join(thisdir, 'test-data/random-20-a.fa')
        outfile = filename + '.out'

        ht.do_threaded_partition(filename)
        subset = ht.do_subset_partition(0, 0)
        n = ht.subset_count_partitions(subset)
        assert n == (0, 99, 0), n

    def test_random_20_b_succ(self):
        ht = khmer.new_hashbits(20, 4**13+1)
        filename = os.path.join(thisdir, 'test-data/random-20-b.fa')
        outfile = filename + '.out'

        ht.do_threaded_partition(filename)
        subset = ht.do_subset_partition(0, 0)
        n = ht.subset_count_partitions(subset)
        assert n == (1, 0, 0), n

    def test_random_20_b_fail(self):
        ht = khmer.new_hashbits(21, 4**13+1)
        filename = os.path.join(thisdir, 'test-data/random-20-b.fa')
        outfile = filename + '.out'

        ht.do_threaded_partition(filename)
        subset = ht.do_subset_partition(0, 0)
        n = ht.subset_count_partitions(subset)
        assert n == (0, 99, 0), n

    def test_random_31_a_succ(self):
        ht = khmer.new_hashbits(31, 4**14+1)
        filename = os.path.join(thisdir, 'test-data/random-31-c.fa')
        outfile = filename + '.out'

        ht.do_threaded_partition(filename)
        subset = ht.do_subset_partition(0, 0)
        n = ht.subset_count_partitions(subset)
        assert n == (1, 0, 0), n

    def test_random_31_b_fail(self):
        ht = khmer.new_hashbits(32, 4**14+1)
        filename = os.path.join(thisdir, 'test-data/random-31-c.fa')
        outfile = filename + '.out'

        ht.do_threaded_partition(filename)
        subset = ht.do_subset_partition(0, 0)
        n = ht.subset_count_partitions(subset)
        assert n == (0, 999, 0), n
