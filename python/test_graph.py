import os
thisdir = os.path.dirname(__file__)
thisdir = os.path.abspath(thisdir)

import khmer

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
        n = ht.do_truncated_partition(filename, outfile, 0)
        assert n == 1, n
        
    def test_simple_30_12(self):
        ht = khmer.new_hashtable(32, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph2.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile, 0)
        assert n == 3, n

class Test_NoConnectMe4(object):
    def test_merge_20_12(self):
        ht = khmer.new_hashtable(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph3.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile, 0)
        assert n == 2, n
        
    def test_merge_32_12(self):
        ht = khmer.new_hashtable(32, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph3.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile, 0)
        assert n == 4, n

class Test_AnotherConnectMe4(object):
    def test_complex_20_12(self):
        ht = khmer.new_hashtable(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph4.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile, 0)
        assert n == 2, n

    def test_complex_31_12(self):
        ht = khmer.new_hashtable(31, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph4.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile, 0)
        assert n == 4, n

    def test_complex_32_12(self):
        ht = khmer.new_hashtable(32, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph4.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile, 0)
        assert n == 5, n


class Test_MoreConnectMe4(object):
    def test_complex5_20_12(self):
        ht = khmer.new_hashtable(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph5.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile, 0)
        assert n == 1, n

    def test_complex5_24_12(self):
        ht = khmer.new_hashtable(30, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph5.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile, 0)
        assert n == 6, n

    def test_complex6_32_12(self):
        ht = khmer.new_hashtable(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph6.fa')
        outfile = os.path.join(thisdir, 'test-trunc.out')
        n = ht.do_truncated_partition(filename, outfile, 0)
        assert n == 103, n

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
        assert pid == 1, pid
        
        print ht.consume(b)
        ppi = ht.find_all_tags(b[:20])
        pid = ht.assign_partition_id(ppi)
        assert pid == 2, pid
        
        print ht.consume(c)
        ppi = ht.find_all_tags(c[:20])
        pid = ht.assign_partition_id(ppi)
        assert pid == 1, pid
