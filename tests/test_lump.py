#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=missing-docstring
import khmer
import screed

import khmer_tst_utils as utils
from nose.plugins.attrib import attr

# Below, 'fakelump.fa' is an artificial data set of 3x1 kb sequences in
# which the last 79 bases are common between the 3 sequences.


@attr('highmem')
def test_fakelump_together():
    fakelump_fa = utils.get_test_data('fakelump.fa')

    ht = khmer.new_hashbits(32, 1e5, 4)
    ht.consume_fasta_and_tag(fakelump_fa)

    subset = ht.do_subset_partition(0, 0)
    ht.merge_subset(subset)

    (n_partitions, n_singletons) = ht.count_partitions()
    assert n_partitions == 1, n_partitions

# try loading stop tags from previously saved


@attr('highmem')
def test_fakelump_stop():
    fakelump_fa = utils.get_test_data('fakelump.fa')
    fakelump_stoptags_txt = utils.get_test_data('fakelump.fa.stoptags.txt')

    ht = khmer.new_hashbits(32, 1e5, 4)
    ht.consume_fasta_and_tag(fakelump_fa)

    for line in open(fakelump_stoptags_txt):
        ht.add_stop_tag(line.strip())

    subset = ht.do_subset_partition(0, 0, True)
    ht.merge_subset(subset)

    (n_partitions, n_singletons) = ht.count_partitions()
    assert n_partitions == 3, n_partitions

# check specific insertion of stop tag


@attr('highmem')
def test_fakelump_stop2():
    fakelump_fa = utils.get_test_data('fakelump.fa')

    ht = khmer.new_hashbits(32, 1e5, 4)
    ht.consume_fasta_and_tag(fakelump_fa)

    ht.add_stop_tag('GGGGAGGGGTGCAGTTGTGACTTGCTCGAGAG')

    subset = ht.do_subset_partition(0, 0, True)
    ht.merge_subset(subset)

    (n_partitions, n_singletons) = ht.count_partitions()
    assert n_partitions == 3, n_partitions

# try repartitioning


@attr('highmem')
def test_fakelump_repartitioning():
    fakelump_fa = utils.get_test_data('fakelump.fa')
    fakelump_fa_foo = utils.get_temp_filename('fakelump.fa.stopfoo')

    ht = khmer.new_hashbits(32, 1e5, 4)
    ht.consume_fasta_and_tag(fakelump_fa)

    subset = ht.do_subset_partition(0, 0)
    ht.merge_subset(subset)

    (n_partitions, n_singletons) = ht.count_partitions()
    assert n_partitions == 1, n_partitions

    # now, break partitions on any k-mer that you see more than once
    # on big excursions, where big excursions are excursions 40 out
    # that encounter more than 82 k-mers.  This should specifically
    # identify our connected sequences in fakelump...

    EXCURSION_DISTANCE = 40
    EXCURSION_KMER_THRESHOLD = 82
    EXCURSION_KMER_COUNT_THRESHOLD = 1
    counting = khmer.new_counting_hash(32, 1e4, 4)

    ht.repartition_largest_partition(None, counting,
                                     EXCURSION_DISTANCE,
                                     EXCURSION_KMER_THRESHOLD,
                                     EXCURSION_KMER_COUNT_THRESHOLD)

    ht.save_stop_tags(fakelump_fa_foo)

    # ok, now re-do everything with these stop tags, specifically.

    ht = khmer.new_hashbits(32, 1e5, 4)
    ht.consume_fasta_and_tag(fakelump_fa)
    ht.load_stop_tags(fakelump_fa_foo)

    subset = ht.do_subset_partition(0, 0, True)
    ht.merge_subset(subset)

    (n_partitions, n_singletons) = ht.count_partitions()
    assert n_partitions == 3, n_partitions


@attr('highmem')
def test_fakelump_load_stop_tags_trunc():
    fakelump_fa = utils.get_test_data('fakelump.fa')
    fakelump_fa_foo = utils.get_temp_filename('fakelump.fa.stopfoo')

    ht = khmer.new_hashbits(32, 1e5, 4)
    ht.consume_fasta_and_tag(fakelump_fa)

    subset = ht.do_subset_partition(0, 0)
    ht.merge_subset(subset)

    (n_partitions, n_singletons) = ht.count_partitions()
    assert n_partitions == 1, n_partitions

    # now, break partitions on any k-mer that you see more than once
    # on big excursions, where big excursions are excursions 40 out
    # that encounter more than 82 k-mers.  This should specifically
    # identify our connected sequences in fakelump...

    EXCURSION_DISTANCE = 40
    EXCURSION_KMER_THRESHOLD = 82
    EXCURSION_KMER_COUNT_THRESHOLD = 1
    counting = khmer.new_counting_hash(32, 4, 4)

    ht.repartition_largest_partition(None, counting,
                                     EXCURSION_DISTANCE,
                                     EXCURSION_KMER_THRESHOLD,
                                     EXCURSION_KMER_COUNT_THRESHOLD)

    ht.save_stop_tags(fakelump_fa_foo)
    data = open(fakelump_fa_foo).read()

    fp = open(fakelump_fa_foo, 'wb')
    fp.write(data[:10])
    fp.close()

    # ok, now try loading these stop tags; should fail.
    ht = khmer.new_hashbits(32, 4, 4)
    ht.consume_fasta_and_tag(fakelump_fa)

    try:
        ht.load_stop_tags(fakelump_fa_foo)
        assert 0, "this test should fail"
    except IOError:
        pass


@attr('highmem')
def test_fakelump_load_stop_tags_notexist():
    fakelump_fa_foo = utils.get_temp_filename('fakelump.fa.stopfoo')

    # ok, now try loading these stop tags; should fail.
    ht = khmer.new_hashbits(32, 4, 4)

    try:
        ht.load_stop_tags(fakelump_fa_foo)
        assert 0, "this test should fail"
    except IOError:
        pass
