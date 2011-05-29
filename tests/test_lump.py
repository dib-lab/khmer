import os
thisdir = os.path.dirname(__file__)
thisdir = os.path.abspath(thisdir)

import khmer
import screed

## Below, 'fakelump.fa' is an artificial data set of 3x1 kb sequences in
## which the last 79 bases are common between the 3 sequences.

def test_fakelump_together():
    fakelump_fa = os.path.join(thisdir, 'test-data/fakelump.fa')

    ht = khmer.new_hashbits(32, 1e7, 4)
    ht.consume_fasta_and_tag(fakelump_fa)

    subset = ht.do_subset_partition(0, 0)
    ht.merge_subset(subset)
    
    (n_partitions, n_singletons) = ht.count_partitions()
    assert n_partitions == 1, n_partitions

# try loading stop tags from previously saved
def test_fakelump_stop():
    fakelump_fa = os.path.join(thisdir, 'test-data/fakelump.fa')
    fakelump_fa_stop = os.path.join(thisdir, 'test-data/fakelump.fa.stoptags')

    ht = khmer.new_hashbits(32, 1e7, 4)
    ht.consume_fasta_and_tag(fakelump_fa)

    ht.load_stop_tags(fakelump_fa_stop)

    subset = ht.do_subset_partition(0, 0, True)
    ht.merge_subset(subset)
    
    (n_partitions, n_singletons) = ht.count_partitions()
    assert n_partitions == 3, n_partitions

# check specific insertion of stop tag
def test_fakelump_stop2():
    fakelump_fa = os.path.join(thisdir, 'test-data/fakelump.fa')

    ht = khmer.new_hashbits(32, 1e7, 4)
    ht.consume_fasta_and_tag(fakelump_fa)

    ht.add_stop_tag('GGGGAGGGGTGCAGTTGTGACTTGCTCGAGAG')

    subset = ht.do_subset_partition(0, 0, True)
    ht.merge_subset(subset)
    
    (n_partitions, n_singletons) = ht.count_partitions()
    assert n_partitions == 3, n_partitions

