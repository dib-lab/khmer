import os
thisdir = os.path.dirname(__file__)
thisdir = os.path.abspath(thisdir)

import khmer

class Test_RandomData(object):
    def test_3_merge_013(self):
        ht = khmer.new_hashtable(20, 4**14+1)
        filename = os.path.join(thisdir, 'test-graph2.fa')

        (total_reads, total_kmers) = ht.consume_fasta_and_tag(filename)
        assert total_reads == 3, total_reads
        
        n_partitions = ht.output_partitions(filename, filename + '.out')
        assert n_partitions == 3, n_partitions        # all singular
        
        x = ht.do_subset_partition(filename, 0, 1)
        ht.merge_subset(x)
        n_partitions = ht.output_partitions(filename, filename + '.out')
        assert n_partitions == 3, n_partitions        # all singular
        
        y = ht.do_subset_partition(filename, 1, 3)
        ht.merge_subset(y)
        
        n_partitions = ht.output_partitions(filename, filename + '.out')
        assert n_partitions == 1, n_partitions        # combined.
    
    def test_3_merge_023(self):
        ht = khmer.new_hashtable(20, 4**14+1)
        filename = os.path.join(thisdir, 'test-graph2.fa')

        (total_reads, total_kmers) = ht.consume_fasta_and_tag(filename)
        assert total_reads == 3, total_reads
        
        n_partitions = ht.output_partitions(filename, filename + '.out')
        assert n_partitions == 3, n_partitions        # all singular
        
        x = ht.do_subset_partition(filename, 0, 2)
        ht.merge_subset(x)
        n_partitions = ht.output_partitions(filename, filename + '.out')
        assert n_partitions == 2, n_partitions        # all singular

        y = ht.do_subset_partition(filename, 2, 3)
        ht.merge_subset(y)
        
        n_partitions = ht.output_partitions(filename, filename + '.out')
        assert n_partitions == 1, n_partitions        # combined.
    
    def test_5_merge_046(self):
        ht = khmer.new_hashtable(20, 4**14+1)
        filename = os.path.join(thisdir, 'test-graph5.fa')

        (total_reads, total_kmers) = ht.consume_fasta_and_tag(filename)
        assert total_reads == 6, total_reads
        
        n_partitions = ht.output_partitions(filename, filename + '.out')
        assert n_partitions == total_reads, n_partitions # all singular
        
        x = ht.do_subset_partition(filename, 0, 4)
        ht.merge_subset(x)
        n_partitions = ht.output_partitions(filename, filename + '.out')
        assert n_partitions == 2, n_partitions

        y = ht.do_subset_partition(filename, 4, 6)
        ht.merge_subset(y)
        
        n_partitions = ht.output_partitions(filename, filename + '.out')
        assert n_partitions == 1, n_partitions        # combined.
    
    def test_random_20_a_succ(self):
        ht = khmer.new_hashtable(20, 4**13+1)
        filename = os.path.join(thisdir, 'test-data/random-20-a.fa')
        outfile = filename + '.out'

        total_reads, _ = ht.consume_fasta_and_tag(filename)
        x = ht.do_subset_partition(filename, 0, total_reads/2)
        ht.merge_subset(x)
        y = ht.do_subset_partition(filename, total_reads/2, total_reads)
        ht.merge_subset(y)

        n_partitions = ht.output_partitions(filename, outfile)

    def test_random_20_a_succ_II(self):
        ht = khmer.new_hashtable(20, 4**13+1)
        filename = os.path.join(thisdir, 'test-data/random-20-a.fa')
        outfile = filename + '.out'

        total_reads, _ = ht.consume_fasta_and_tag(filename)
        x = ht.do_subset_partition(filename, 0, total_reads/2)
        y = ht.do_subset_partition(filename, total_reads/2, total_reads)
        ht.merge_subset(x)
        ht.merge_subset(y)

        n_partitions = ht.output_partitions(filename, outfile)

    def test_random_20_a_succ_III(self):
        ht = khmer.new_hashtable(20, 4**13+1)
        filename = os.path.join(thisdir, 'test-data/random-20-a.fa')
        outfile = filename + '.out'

        total_reads, _ = ht.consume_fasta_and_tag(filename)
        x = ht.do_subset_partition(filename, 0, total_reads/2)
        y = ht.do_subset_partition(filename, total_reads/2, total_reads)
        ht.merge_subset(y)
        ht.merge_subset(x)

        n_partitions = ht.output_partitions(filename, outfile)

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
