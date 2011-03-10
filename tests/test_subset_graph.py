import os
thisdir = os.path.dirname(__file__)
thisdir = os.path.abspath(thisdir)

import khmer
import screed

class Test_RandomData(object):
    def test_3_merge_013(self):
        ht = khmer.new_hashbits(20, 4**14+1)
        filename = os.path.join(thisdir, 'test-graph2.fa')

        (total_reads, total_kmers) = ht.consume_fasta_and_tag(filename)
        assert total_reads == 3, total_reads
        
        (a, b, c) = ht.divide_tags_into_subsets(1)
        
        x = ht.do_subset_partition(a, a)
        ht.merge_subset(x)
        
        y = ht.do_subset_partition(b, 0)
        ht.merge_subset(y)
        
        n_partitions = ht.output_partitions(filename, filename + '.out')
        assert n_partitions == 1, n_partitions        # combined.
    
    def test_3_merge_023(self):
        ht = khmer.new_hashbits(20, 4**14+1)
        filename = os.path.join(thisdir, 'test-graph2.fa')

        (total_reads, total_kmers) = ht.consume_fasta_and_tag(filename)
        assert total_reads == 3, total_reads
        
        (a, b, c) = ht.divide_tags_into_subsets(1)
        
        x = ht.do_subset_partition(b, c)
        ht.merge_subset(x)

        y = ht.do_subset_partition(a, b)
        ht.merge_subset(y)
        
        n_partitions = ht.output_partitions(filename, filename + '.out')
        assert n_partitions == 1, n_partitions        # combined.
    
    def test_5_merge_046(self):
        ht = khmer.new_hashbits(20, 4**14+1)
        filename = os.path.join(thisdir, 'test-graph5.fa')

        (total_reads, total_kmers) = ht.consume_fasta_and_tag(filename)
        assert total_reads == 6, total_reads
        
        divvy = ht.divide_tags_into_subsets(1)
        
        x = ht.do_subset_partition(divvy[0], divvy[4])
        ht.merge_subset(x)

        y = ht.do_subset_partition(divvy[4], 0)
        ht.merge_subset(y)
        
        n_partitions = ht.output_partitions(filename, filename + '.out')
        assert n_partitions == 1, n_partitions        # combined.
    
    def test_random_20_a_succ(self):
        ht = khmer.new_hashbits(20, 4**13+1)
        filename = os.path.join(thisdir, 'test-data/random-20-a.fa')
        outfile = filename + '.out'

        total_reads, _ = ht.consume_fasta_and_tag(filename)

        subset_size = total_reads / 2 + total_reads % 2;
        divvy = ht.divide_tags_into_subsets(subset_size)
        assert len(divvy) == 4

        x = ht.do_subset_partition(divvy[0], divvy[2])
        ht.merge_subset(x)
        y = ht.do_subset_partition(divvy[2], 0)
        ht.merge_subset(y)

        n_partitions = ht.output_partitions(filename, outfile)
        assert n_partitions == 1, n_partitions

    def test_random_20_a_succ_II(self):
        ht = khmer.new_hashbits(20, 4**13+1)
        filename = os.path.join(thisdir, 'test-data/random-20-a.fa')
        outfile = filename + '.out'

        total_reads, _ = ht.consume_fasta_and_tag(filename)

        subset_size = total_reads / 2 + total_reads % 2;
        divvy = ht.divide_tags_into_subsets(subset_size)
        assert len(divvy) == 4

        x = ht.do_subset_partition(divvy[0], divvy[2])
        y = ht.do_subset_partition(divvy[2], 0)
        ht.merge_subset(x)
        ht.merge_subset(y)

        n_partitions = ht.output_partitions(filename, outfile)
        assert n_partitions == 1, n_partitions

    def test_random_20_a_succ_III(self):
        ht = khmer.new_hashbits(20, 4**13+1)
        filename = os.path.join(thisdir, 'test-data/random-20-a.fa')
        outfile = filename + '.out'

        total_reads, _ = ht.consume_fasta_and_tag(filename)

        subset_size = total_reads / 2 + total_reads % 2;
        divvy = ht.divide_tags_into_subsets(subset_size)
        assert len(divvy) == 4, len(divvy)

        x = ht.do_subset_partition(divvy[0], divvy[2])
        y = ht.do_subset_partition(divvy[2], 0)

        ht._validate_subset_partitionmap(x)
        ht._validate_subset_partitionmap(y)
        
        ht.merge_subset(y)
        ht.merge_subset(x)

        n_partitions = ht.output_partitions(filename, outfile)
        assert n_partitions == 1, n_partitions

    def test_random_20_a_succ_IV(self):
        ht = khmer.new_hashbits(20, 4**13+1)
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
        ht = khmer.new_hashbits(20, 4**13+1)
        filename = os.path.join(thisdir, 'test-data/random-20-a.fa')
        savefile_ht = filename + '.ht'
        savefile_tags = filename + '.tags'
        outfile = filename + '.out'

        total_reads, _ = ht.consume_fasta_and_tag(filename)

        ht.save(savefile_ht);
        ht.save_tagset(savefile_tags);

        del ht
        ht = khmer.new_hashbits(20, 4**13+1)

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

class Test_SaveLoadPmap(object):
    def test_save_load_merge(self):
        ht = khmer.new_hashbits(20, 4**14+1)
        filename = os.path.join(thisdir, 'test-graph2.fa')

        (total_reads, total_kmers) = ht.consume_fasta_and_tag(filename)
        assert total_reads == 3, total_reads

        divvy = ht.divide_tags_into_subsets(1)
        print divvy
        (a, b, c) = divvy
        
        x = ht.do_subset_partition(a, b)
        ht.save_subset_partitionmap(x, 'x.pmap')
        del x

        y = ht.do_subset_partition(b, 0)
        ht.save_subset_partitionmap(y, 'y.pmap')
        del y

        a = ht.load_subset_partitionmap('x.pmap')
        b = ht.load_subset_partitionmap('y.pmap')

        ht.merge_subset(a)
        ht.merge_subset(b)
        
        n_partitions = ht.output_partitions(filename, filename + '.out')
        assert n_partitions == 1, n_partitions        # combined.
        
    def test_save_load_merge_2(self):
        ht = khmer.new_hashbits(20, 4**14+1)
        filename = os.path.join(thisdir, 'test-data/random-20-a.fa')

        (total_reads, total_kmers) = ht.consume_fasta_and_tag(filename)

        subset_size = total_reads / 2 + total_reads % 2;
        divvy = ht.divide_tags_into_subsets(subset_size)
        
        x = ht.do_subset_partition(divvy[0], divvy[1])
        ht.save_subset_partitionmap(x, 'x.pmap')
        del x

        y = ht.do_subset_partition(divvy[1], 0)
        ht.save_subset_partitionmap(y, 'y.pmap')
        del y

        a = ht.load_subset_partitionmap('x.pmap')
        b = ht.load_subset_partitionmap('y.pmap')

        ht.merge_subset(a)
        ht.merge_subset(b)
        
        n_partitions = ht.output_partitions(filename, filename + '.out')
        assert n_partitions == 1, n_partitions        # combined.

    def test_save_merge_from_disk(self):
        ht = khmer.new_hashbits(20, 4**14+1)
        filename = os.path.join(thisdir, 'test-graph2.fa')

        (total_reads, total_kmers) = ht.consume_fasta_and_tag(filename)
        assert total_reads == 3, total_reads

        divvy = ht.divide_tags_into_subsets(1)
        print divvy
        (a, b, c) = divvy
        
        x = ht.do_subset_partition(a, b)
        ht.save_subset_partitionmap(x, 'x.pmap')
        del x

        y = ht.do_subset_partition(b, 0)
        ht.save_subset_partitionmap(y, 'y.pmap')
        del y

        ht.merge_subset_from_disk('x.pmap')
        ht.merge_subset_from_disk('y.pmap')
        
        n_partitions = ht.output_partitions(filename, filename + '.out')
        assert n_partitions == 1, n_partitions        # combined.
        
    def test_save_merge_from_disk_2(self):
        ht = khmer.new_hashbits(20, 4**14+1)
        filename = os.path.join(thisdir, 'test-data/random-20-a.fa')

        (total_reads, total_kmers) = ht.consume_fasta_and_tag(filename)

        subset_size = total_reads / 2 + total_reads % 2;
        divvy = ht.divide_tags_into_subsets(subset_size)
        
        x = ht.do_subset_partition(divvy[0], divvy[1])
        ht.save_subset_partitionmap(x, 'x.pmap')
        del x

        y = ht.do_subset_partition(divvy[1], 0)
        ht.save_subset_partitionmap(y, 'y.pmap')
        del y

        ht.merge_subset_from_disk('x.pmap')
        ht.merge_subset_from_disk('y.pmap')
        
        n_partitions = ht.output_partitions(filename, filename + '.out')
        assert n_partitions == 1, n_partitions        # combined.

def test_output_partitions():
    filename = os.path.join(thisdir, 'test-data/test-output-partitions.fa')

    ht = khmer.new_hashbits(10, 1, 1)
    ht.set_partition_id('TTAGGACTGC', 2)
    ht.set_partition_id('TGCGTTTCAA', 3)
    ht.set_partition_id('ATACTGTAAA', 4)

    try:
        os.unlink(filename + '.part')
    except OSError:
        pass
    assert not os.path.exists(filename + '.part')
    ht.output_partitions(filename, filename + '.part')

    data = open(filename + '.part').read()
    assert len(data)

    records = [ r for r in screed.open(filename + '.part') ]
    names = [ r.name for r in records ]
    parts = [ n.rsplit('\t', 1)[1] for n in names ]

    assert parts[0] == '2'
    assert parts[1] == '3'
    assert parts[2] == '4'

test_output_partitions.runme = True

def test_tiny_real_partitions():
    filename = os.path.join(thisdir, 'test-data/real-partition-tiny.fa')
    
    ht = khmer.new_hashbits(32, 8e7, 4)
    ht.consume_fasta_and_tag(filename)
    
    subset = ht.do_subset_partition(0, 0)
    ht.merge_subset(subset)

    try:
        os.unlink(filename + '.part')
    except OSError:
        pass
    
    assert not os.path.exists(filename + '.part')
    ht.output_partitions(filename, filename + '.part')

    data = open(filename + '.part').read()
    assert len(data)
    
    records = [ r for r in screed.open(filename + '.part') ]
    names = [ r.name for r in records ]
    parts = [ n.rsplit('\t', 1)[1] for n in names ]

    assert len(parts) == 2, len(parts)
    assert len(set(parts)) == 1
    assert set(parts) != set(['0'])

test_tiny_real_partitions.runme = True

def test_small_real_partitions():
    filename = os.path.join(thisdir, 'test-data/real-partition-small.fa')
    
    ht = khmer.new_hashbits(32, 8e7, 4)
    ht.consume_fasta_and_tag(filename)

    subset = ht.do_subset_partition(0, 0)
    ht.merge_subset(subset)

    try:
        os.unlink(filename + '.part')
    except OSError:
        pass
    
    assert not os.path.exists(filename + '.part')
    ht.output_partitions(filename, filename + '.part')

    data = open(filename + '.part').read()
    assert len(data)
    
    records = [ r for r in screed.open(filename + '.part') ]
    names = [ r.name for r in records ]
    parts = [ n.rsplit('\t', 1)[1] for n in names ]

    assert len(parts) == 6, len(parts)
    assert len(set(parts)) == 1
    assert set(parts) != set(['0'])

test_small_real_partitions.runme = True
