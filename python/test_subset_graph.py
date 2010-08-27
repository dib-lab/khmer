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

        (a, b, c) = ht.divide_tags_into_subsets(1)
        
        x = ht.do_subset_partition(a, a)
        ht.merge_subset(x)
        n_partitions = ht.output_partitions(filename, filename + '.out')
        assert n_partitions == 3, n_partitions        # all singular
        
        y = ht.do_subset_partition(b, 0)
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

        (a, b, c) = ht.divide_tags_into_subsets(1)
        
        x = ht.do_subset_partition(a, b)
        ht.merge_subset(x)
        n_partitions = ht.output_partitions(filename, filename + '.out')
        assert n_partitions == 2, n_partitions        # all singular

        y = ht.do_subset_partition(b, c)
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

        divvy = ht.divide_tags_into_subsets(1)
        
        x = ht.do_subset_partition(divvy[0], divvy[4])
        ht.merge_subset(x)
        n_partitions = ht.output_partitions(filename, filename + '.out')
        assert n_partitions == 2, n_partitions

        y = ht.do_subset_partition(divvy[4], 0)
        ht.merge_subset(y)
        
        n_partitions = ht.output_partitions(filename, filename + '.out')
        assert n_partitions == 1, n_partitions        # combined.
    
    def test_random_20_a_succ(self):
        ht = khmer.new_hashtable(20, 4**13+1)
        filename = os.path.join(thisdir, 'test-data/random-20-a.fa')
        outfile = filename + '.out'

        total_reads, _ = ht.consume_fasta_and_tag(filename)

        subset_size = total_reads / 2 + total_reads % 2;
        divvy = ht.divide_tags_into_subsets(subset_size)
        assert len(divvy) == 2

        x = ht.do_subset_partition(divvy[0], divvy[1])
        ht.merge_subset(x)
        y = ht.do_subset_partition(divvy[0], 0)
        ht.merge_subset(y)

        n_partitions = ht.output_partitions(filename, outfile)
        assert n_partitions == 1, n_partitions

    def test_random_20_a_succ_II(self):
        ht = khmer.new_hashtable(20, 4**13+1)
        filename = os.path.join(thisdir, 'test-data/random-20-a.fa')
        outfile = filename + '.out'

        total_reads, _ = ht.consume_fasta_and_tag(filename)

        subset_size = total_reads / 2 + total_reads % 2;
        divvy = ht.divide_tags_into_subsets(subset_size)
        assert len(divvy) == 2

        x = ht.do_subset_partition(divvy[0], divvy[1])
        y = ht.do_subset_partition(divvy[1], 0)
        ht.merge_subset(x)
        ht.merge_subset(y)

        n_partitions = ht.output_partitions(filename, outfile)
        assert n_partitions == 1, n_partitions

    def test_random_20_a_succ_III(self):
        ht = khmer.new_hashtable(20, 4**13+1)
        filename = os.path.join(thisdir, 'test-data/random-20-a.fa')
        outfile = filename + '.out'

        total_reads, _ = ht.consume_fasta_and_tag(filename)

        subset_size = total_reads / 2 + total_reads % 2;
        divvy = ht.divide_tags_into_subsets(subset_size)
        assert len(divvy) == 2

        x = ht.do_subset_partition(divvy[0], divvy[1])
        y = ht.do_subset_partition(divvy[1], 0)

        ht._validate_subset_partitionmap(x)
        ht._validate_subset_partitionmap(y)
        
        ht.merge_subset(y)
        ht.merge_subset(x)

        n_partitions = ht.output_partitions(filename, outfile)
        assert n_partitions == 1, n_partitions

    def test_random_20_a_succ_IV(self):
        ht = khmer.new_hashtable(20, 4**13+1)
        filename = os.path.join(thisdir, 'test-data/random-20-a.fa')
        outfile = filename + '.out'

        total_reads, _ = ht.consume_fasta_and_tag(filename)
        subsets = []

        divvy = ht.divide_tags_into_subsets(1)
        divvy.append(0)
        for i in range(len(divvy) - 1):
            x = ht.do_subset_partition(divvy[i], divvy[i+1])
            subsets.append(x)

        for x in reversed(subsets):
            ht.merge_subset(x)
            
        n_partitions = ht.output_partitions(filename, outfile)
        assert n_partitions == 1, n_partitions
        
    def test_random_20_a_succ_IV_save(self):
        ht = khmer.new_hashtable(20, 4**13+1)
        filename = os.path.join(thisdir, 'test-data/random-20-a.fa')
        savefile_ht = filename + '.ht'
        savefile_tags = filename + '.tags'
        outfile = filename + '.out'

        total_reads, _ = ht.consume_fasta_and_tag(filename)

        ht.save(savefile_ht);
        ht.save_tagset(savefile_tags);

        del ht
        ht = khmer.new_hashtable(20, 4**13+1)

        ht.load(savefile_ht);
        ht.load_tagset(savefile_tags);
        
        divvy = ht.divide_tags_into_subsets(1)
        divvy.append(0)
        
        subsets = []
        for i in range(len(divvy) - 1):
            x = ht.do_subset_partition(divvy[i], divvy[i+1])
            subsets.append(x)

        for x in reversed(subsets):
            ht.merge_subset(x)
            
        n_partitions = ht.output_partitions(filename, outfile)
        assert n_partitions == 1, n_partitions

class Test_Surrendered(object):
    def test_surrendered_subset(self):
        ht = khmer.new_hashtable(32, 4**15+1)

        filename = os.path.join(thisdir, '../data/100k-surrendered.fa')
        total_reads, _ = ht.consume_fasta_and_tag(filename)
        subset = ht.do_subset_partition(0, total_reads)
        ht.merge_subset(subset)

        n_partitions, n_unassigned, n_surrendered = ht.count_partitions()

        assert n_partitions == 16, n_partitions
        assert n_unassigned == 0, n_unassigned
        assert n_surrendered == 16, n_surrendered

    
    def test_surrendered_subset_2(self):
        ht = khmer.new_hashtable(32, 4**15+1)

        filename = os.path.join(thisdir, '../data/100k-surrendered.fa')
        total_reads, _ = ht.consume_fasta_and_tag(filename)
        subset1 = ht.do_subset_partition(0, total_reads / 2)
        subset2 = ht.do_subset_partition(total_reads / 2, total_reads)
        ht.merge_subset(subset1)
        ht.merge_subset(subset2)

        n_partitions, n_unassigned, n_surrendered = ht.count_partitions()

        assert n_partitions == 16, n_partitions
        assert n_unassigned == 0, n_unassigned
        assert n_surrendered == 16, n_surrendered

    def test_save_load_merge(self):
        ht = khmer.new_hashtable(20, 4**14+1)
        filename = os.path.join(thisdir, 'test-graph2.fa')

        (total_reads, total_kmers) = ht.consume_fasta_and_tag(filename)
        assert total_reads == 3, total_reads

        divvy = ht.divide_tags_into_subsets(1)
        (a, b, c) = divvy
        
        n_partitions = ht.output_partitions(filename, filename + '.out')
        assert n_partitions == 3, n_partitions        # all singular
        
        x = ht.do_subset_partition(a, b)
        ht.save_subset_partitionmap(x, 'x.pmap', 'x.surrender')
        del x

        y = ht.do_subset_partition(b, 0)
        ht.save_subset_partitionmap(y, 'y.pmap', 'y.surrender')
        del y

        a = ht.load_subset_partitionmap('x.pmap', 'x.surrender')
        b = ht.load_subset_partitionmap('y.pmap', 'y.surrender')

        ht.merge_subset(a)
        ht.merge_subset(b)
        
        n_partitions = ht.output_partitions(filename, filename + '.out')
        assert n_partitions == 1, n_partitions        # combined.
        
    def test_save_load_merge_2(self):
        ht = khmer.new_hashtable(20, 4**14+1)
        filename = os.path.join(thisdir, 'test-data/random-20-a.fa')

        (total_reads, total_kmers) = ht.consume_fasta_and_tag(filename)

        subset_size = total_reads / 2 + total_reads % 2;
        divvy = ht.divide_tags_into_subsets(subset_size)
        
        x = ht.do_subset_partition(divvy[0], divvy[1])
        ht.save_subset_partitionmap(x, 'x.pmap', 'x.surrender')
        del x

        y = ht.do_subset_partition(divvy[1], 0)
        ht.save_subset_partitionmap(y, 'y.pmap', 'y.surrender')
        del y

        a = ht.load_subset_partitionmap('x.pmap', 'x.surrender')
        b = ht.load_subset_partitionmap('y.pmap', 'y.surrender')

        ht.merge_subset(a)
        ht.merge_subset(b)
        
        n_partitions = ht.output_partitions(filename, filename + '.out')
        assert n_partitions == 1, n_partitions        # combined.

    def test_surrendered_save_load(self):
        ht = khmer.new_hashtable(32, 4**15+1)

        filename = os.path.join(thisdir, '../data/100k-surrendered.fa')
        total_reads, _ = ht.consume_fasta_and_tag(filename)
        subset = ht.do_subset_partition(0, 0)

        ht.save_subset_partitionmap(subset, 'aaa.pmap', 'aaa.surr')
        del subset
        
        subset = ht.load_subset_partitionmap('aaa.pmap', 'aaa.surr')
        ht.merge_subset(subset)

        n_partitions, n_unassigned, n_surrendered = ht.count_partitions()

        assert n_partitions == 16, n_partitions
        assert n_unassigned == 0, n_unassigned
        assert n_surrendered == 16, n_surrendered

        ht.output_partitions(filename, filename + '.out')

        for line in open(filename + '.out'):
            line = line.strip()
            if line.startswith('>'):
                assert line.endswith('*'), line

