import os
thisdir = os.path.dirname(__file__)
thisdir = os.path.abspath(thisdir)

import khmer

class Test_SimpleConnectMe(object):
    def test_simple_20_12(self):
        ht = khmer.new_hashtable(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph2.fa')
        n = ht.do_partition(filename, "")
        assert n == 1, x
        
    def test_simple_30_12(self):
        ht = khmer.new_hashtable(32, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph2.fa')
        n = ht.do_partition(filename, "")
        assert n == 3, x

class Test_NoConnectMe(object):
    def test_merge_20_12(self):
        ht = khmer.new_hashtable(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph3.fa')
        n = ht.do_partition(filename, "")
        assert n == 2, x
        
    def test_merge_32_12(self):
        ht = khmer.new_hashtable(32, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph3.fa')
        n = ht.do_partition(filename, "")
        assert n == 4, x

class Test_AnotherConnectMe(object):
    def test_complex_20_12(self):
        ht = khmer.new_hashtable(20, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph4.fa')
        n = ht.do_partition(filename, "")
        assert n == 2, x

    def test_complex_31_12(self):
        ht = khmer.new_hashtable(31, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph4.fa')
        n = ht.do_partition(filename, "")
        assert n == 4, x

    def test_complex_32_12(self):
        ht = khmer.new_hashtable(32, 4**12+1)

        filename = os.path.join(thisdir, 'test-graph4.fa')
        n = ht.do_partition(filename, "")
        assert n == 5, x
