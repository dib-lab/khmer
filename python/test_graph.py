import os
thisdir = os.path.dirname(__file__)
thisdir = os.path.abspath(thisdir)

import khmer

class Test_PartitionCount(object):
    def test_simple_20_12(self):
        ht = khmer.new_hashtable(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph2.fa')
        ht.do_truncated_partition(filename, filename + '.out')
        n_partitions, n_unassigned, n_surrendered = ht.count_partitions()

        assert n_partitions == 1
        assert n_unassigned == 0
        assert n_surrendered == 0
        
    def test_simple_30_12(self):
        ht = khmer.new_hashtable(32, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph2.fa')
        ht.do_truncated_partition(filename, filename + '.out')
        n_partitions, n_unassigned, n_surrendered = ht.count_partitions()

        assert n_partitions == 0
        assert n_unassigned == 3
        assert n_surrendered == 0

    def test_surrendered(self):
        ht = khmer.new_hashtable(32, 4**15+1)

        filename = os.path.join(thisdir, '../data/100k-surrendered.fa')
        ht.do_truncated_partition(filename, filename + '.out')
        n_partitions, n_unassigned, n_surrendered = ht.count_partitions()

        assert n_partitions == 16, n_partitions
        assert n_unassigned == 0, n_unassigned
        assert n_surrendered == 16, n_surrendered

class Test_SimpleConnectMe(object):
    def test_simple_20_12(self):
        ht = khmer.new_hashtable(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph2.fa')
        n = ht.do_exact_partition(filename)
        assert n == 1, n
        
    def test_simple_30_12(self):
        ht = khmer.new_hashtable(32, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph2.fa')
        n = ht.do_exact_partition(filename)
        assert n == 3, n

class Test_NoConnectMe(object):
    def test_merge_20_12(self):
        ht = khmer.new_hashtable(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph3.fa')
        n = ht.do_exact_partition(filename)
        assert n == 2, n
        
    def test_merge_32_12(self):
        ht = khmer.new_hashtable(32, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph3.fa')
        n = ht.do_exact_partition(filename)
        assert n == 4, n

class Test_AnotherConnectMe(object):
    def test_complex_20_12(self):
        ht = khmer.new_hashtable(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph4.fa')
        n = ht.do_exact_partition(filename)
        assert n == 2, n

    def test_complex_31_12(self):
        ht = khmer.new_hashtable(31, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph4.fa')
        n = ht.do_exact_partition(filename)
        assert n == 4, n

    def test_complex_32_12(self):
        ht = khmer.new_hashtable(32, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph4.fa')
        n = ht.do_exact_partition(filename)
        assert n == 5, n

class Test_MoreConnectMe(object):
    def test_complex5_20_12(self):
        ht = khmer.new_hashtable(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph5.fa')
        n = ht.do_exact_partition(filename)
        assert n == 1, n

    def test_complex5_24_12(self):
        ht = khmer.new_hashtable(30, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph5.fa')
        n = ht.do_exact_partition(filename)
        assert n == 6, n

    def test_complex6_32_12(self):
        ht = khmer.new_hashtable(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph6.fa')
        n = ht.do_exact_partition(filename)
        assert n == 103, n

### do_truncated_partition

class Test_SimpleConnectMe4(object):
    def test_simple_20_12(self):
        ht = khmer.new_hashtable(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph2.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out') # @CTB use tempfile
        n = ht.do_truncated_partition(filename, outfile)
        assert n == 1, n
        
    def test_simple_30_12(self):
        ht = khmer.new_hashtable(32, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph2.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile)
        assert n == 3, n

class Test_NoConnectMe4(object):
    def test_merge_20_12(self):
        ht = khmer.new_hashtable(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph3.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile)
        assert n == 2, n
        
    def test_merge_32_12(self):
        ht = khmer.new_hashtable(32, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph3.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile)
        assert n == 4, n

class Test_AnotherConnectMe4(object):
    def test_complex_20_12(self):
        ht = khmer.new_hashtable(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph4.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile)
        assert n == 2, n

    def test_complex_31_12(self):
        ht = khmer.new_hashtable(31, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph4.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile)
        assert n == 4, n

    def test_complex_32_12(self):
        ht = khmer.new_hashtable(32, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph4.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile)
        assert n == 5, n


class Test_MoreConnectMe4(object):
    def test_complex5_20_12(self):
        ht = khmer.new_hashtable(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph5.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile)
        assert n == 1, n

    def test_complex5_24_12(self):
        ht = khmer.new_hashtable(30, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph5.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile)
        assert n == 6, n

    def test_complex6_32_12(self):
        ht = khmer.new_hashtable(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph6.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile)
        assert n == 103, n

    def test_complex6_32_12_save(self):
        ht = khmer.new_hashtable(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph6.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile)
        assert n == 103, n
        ht._validate_partitionmap()

        o1 = os.path.join(thisdir, 'xx.pmap')
        o2 = os.path.join(thisdir, 'xx.surrender')
        ht.save_partitionmap(o1, o2)
        ht.load_partitionmap(o1, o2)

        n = ht.do_truncated_partition(filename, outfile)
        assert n == 103, n
        ht._validate_partitionmap()

###

class Test_PythonAPI(object):
    def test_ordered_connect(self):
        ht = khmer.new_hashtable(20, 4**15+1)

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
        assert pid == 1, pid

###

class Test_RandomData(object):
    def test_random_20_a_succ(self):
        ht = khmer.new_hashtable(20, 4**13+1)
        filename = os.path.join(thisdir, 'test-data/random-20-a.fa')
        outfile = filename + '.out'

        n = ht.do_truncated_partition(filename, outfile)
        assert n == 1, n

    def test_random_20_a_fail(self):
        ht = khmer.new_hashtable(21, 4**13+1)
        filename = os.path.join(thisdir, 'test-data/random-20-a.fa')
        outfile = filename + '.out'

        n = ht.do_truncated_partition(filename, outfile)
        assert n == 99, n

    def test_random_20_b_succ(self):
        ht = khmer.new_hashtable(20, 4**13+1)
        filename = os.path.join(thisdir, 'test-data/random-20-b.fa')
        outfile = filename + '.out'

        n = ht.do_truncated_partition(filename, outfile)
        assert n == 1, n

    def test_random_20_b_fail(self):
        ht = khmer.new_hashtable(21, 4**13+1)
        filename = os.path.join(thisdir, 'test-data/random-20-b.fa')
        outfile = filename + '.out'

        n = ht.do_truncated_partition(filename, outfile)
        assert n == 99, n

    def test_random_31_a_succ(self):
        ht = khmer.new_hashtable(31, 4**14+1)
        filename = os.path.join(thisdir, 'test-data/random-31-c.fa')
        outfile = filename + '.out'

        n = ht.do_truncated_partition(filename, outfile)
        assert n == 1, n

    def test_random_31_b_fail(self):
        ht = khmer.new_hashtable(32, 4**14+1)
        filename = os.path.join(thisdir, 'test-data/random-31-c.fa')
        outfile = filename + '.out'

        n = ht.do_truncated_partition(filename, outfile)
        assert n == 999, n
