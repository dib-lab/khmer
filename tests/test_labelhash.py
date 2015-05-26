#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# pylint: disable=missing-docstring,protected-access
import khmer
from khmer import LabelHash
from screed.fasta import fasta_iter
import screed

import khmer_tst_utils as utils
from nose.plugins.attrib import attr


def teardown():
    utils.cleanup()

#
# @camillescott TODO: more tests!
#  * thread-safety


@attr('linux')
def test_toobig():
    try:
        lh = LabelHash(20, 1e13, 1)
        assert 0, "This should fail."
    except MemoryError as err:
        print str(err)


def test_error_create():
    from khmer import _LabelHash
    try:
        lh = _LabelHash(None)
        assert 0, "This should fail."
    except ValueError as err:
        print str(err)


def test_n_labels():
    lh = LabelHash(20, 1e7, 4)
    filename = utils.get_test_data('test-labels.fa')
    lh.consume_fasta_and_tag_with_labels(filename)

    print lh.n_labels()
    assert lh.n_labels() == 4


def test_get_label_dict():
    lb = LabelHash(20, 1e7, 4)
    filename = utils.get_test_data('test-labels.fa')
    lb.consume_fasta_and_tag_with_labels(filename)

    labels = lb.get_label_dict()
    expected = [0, 1, 2, 3]
    for e_label in expected:
        assert e_label in labels
    for a_label in labels:
        assert a_label in expected


def test_get_label_dict_save_load():
    lb_pre = LabelHash(20, 1e7, 4)
    filename = utils.get_test_data('test-labels.fa')
    lb_pre.consume_fasta_and_tag_with_labels(filename)

    # save labels to a file
    savepath = utils.get_temp_filename('saved.labels')
    lb_pre.save_labels_and_tags(savepath)

    # trash the old LabelHash
    del lb_pre

    # create new, load labels & tags
    lb = LabelHash(20, 1e7, 4)
    lb.load_labels_and_tags(savepath)

    labels = lb.get_label_dict()
    expected = [0, 1, 2, 3]
    for e_label in expected:
        assert e_label in labels
    for a_label in labels:
        assert a_label in expected


def test_get_tag_labels():
    lb = LabelHash(20, 1e7, 4)
    filename = utils.get_test_data('single-read.fq')
    lb.consume_fasta_and_tag_with_labels(filename)
    tag = 173473779682

    labels = lb.get_tag_labels(tag)
    assert len(labels) == 1
    assert labels.pop() == 0


def test_consume_fasta_and_tag_with_labels():
    lb = LabelHash(20, 1e7, 4)
    read_1 = 'ACGTAACCGGTTAAACCCGGGTTTAAAACCCCGGGGTTTT'
    filename = utils.get_test_data('test-transcript.fa')

    total_reads, n_consumed = lb.consume_fasta_and_tag_with_labels(filename)
    print "doing get"
    assert lb.graph.get(read_1[:20])
    assert total_reads == 3
    print "doing n_labels"
    print lb.n_labels()
    print "doing label dict"
    print lb.get_label_dict()
    print "get tagset"
    for tag in lb.graph.get_tagset():
        print "forward hash"
        print tag, khmer.forward_hash(tag, 20)
    for record in screed.open(filename):
        print "Sweeping tags"
        print lb.sweep_tag_neighborhood(record.sequence, 40)
        print "Sweeping labels..."
        print lb.sweep_label_neighborhood(record.sequence, 40)
    assert lb.n_labels() == 3


def test_consume_partitioned_fasta_and_tag_with_labels():
    lb = LabelHash(20, 1e7, 4)
    filename = utils.get_test_data('real-partition-small.fa')

    total_reads, n_consumed = lb.consume_partitioned_fasta_and_tag_with_labels(
        filename)
    labels = set()
    for record in screed.open(filename):
        seq = record.sequence
        labels.update(lb.sweep_label_neighborhood(seq, 0, False, False))
    # print lb.n_labels()
    # print labels
    assert len(labels) == 1
    assert labels.pop() == 2
    assert lb.n_labels() == 1


def test_consume_sequence_and_tag_with_labels():
    lb = LabelHash(20, 1e6, 4)
    label = 0
    sequence = 'ATGCATCGATCGATCGATCGATCGATCGATCGATCGATCG'

    n_consumed = lb.consume_sequence_and_tag_with_labels(sequence, label)
    labels = set()
    labels.update(lb.sweep_label_neighborhood(sequence))

    assert label in labels
    assert len(labels) == 1


def test_sweep_tag_neighborhood():
    lb = LabelHash(20, 1e7, 4)
    filename = utils.get_test_data('single-read.fq')
    lb.graph.consume_fasta_and_tag(filename)

    tags = lb.sweep_tag_neighborhood('CAGGCGCCCACCACCGTGCCCTCCAACCTGATGGT')
    assert len(tags) == 1
    assert tags.pop() == 173473779682


def test_sweep_label_neighborhood():
    lb = LabelHash(20, 1e7, 4)
    filename = utils.get_test_data('single-read.fq')
    lb.consume_fasta_and_tag_with_labels(filename)

    labels = lb.sweep_label_neighborhood('CAGGCGCCCACCACCGTGCCCTCCAACCTGATGGT')
    assert len(labels) == 1
    assert labels.pop() == 0

'''
* The test data set as four reads: A, B, C, and D
* Overlaps are A <-> B <-> C, with D on its own
* Thus, traversing from A should find labels from A and B,
  traversing from B should find labels from A, B, and C,
  and traversing from C should find labels from B and C
'''


def test_label_tag_correctness():
    lb = LabelHash(20, 1e7, 4)
    filename = utils.get_test_data('test-labels.fa')
    lb.consume_fasta_and_tag_with_labels(filename)

    # read A
    labels = lb.sweep_label_neighborhood(
        'ATCGTGTAAGCTATCGTAATCGTAAGCTCTGCCTAGAGCTAGGCTAGGCTCTGCCTAGAG'
        'CTAGGCTAGGTGTGCTCTGCCTAGAGCTAGGCTAGGTGT')
    print lb.sweep_tag_neighborhood(
        'TTCGTGTAAGCTATCGTAATCGTAAGCTCTGCCTAGAGCTAGGCTAGGCTCTGCCTAGAG'
        'CTAGGCTAGGTGTGCTCTGCTAGAGCTAGGCTAGGTGT')
    print labels
    print len('ATCGTGTAAGCTATCGTAATCGTAAGCTCTGCCTAGAGCTAGGCTAG') - 19
    assert len(labels) == 2
    assert 0 in labels
    assert 1 in labels

    # read B
    labels = lb.sweep_label_neighborhood(
        'GCGTAATCGTAAGCTCTGCCTAGAGCTAGGCTAGCTCTGCCTAGAGCTAGGCTAGGTGTTGGGGATAG'
        'ATAGATAGATGACCTAGAGCTAGGCTAGGTGTTGGGGATAGATAGATAGATGA')
    print labels
    assert len(labels) == 3
    assert 0 in labels
    assert 1 in labels
    assert 2 in labels

    # read C
    labels = lb.sweep_label_neighborhood(
        'TGGGATAGATAGATAGATGACCTAGAGCTAGGCTAGGTGTTGGGGATAGATAGATAGATGACCTAGAG'
        'CTAGGCTAGGTGTTGGGGATAGATAGATAGATGAGTTGGGGATAGATAGATAGATGAGTGTAGATCCA'
        'ACAACACATACA')
    print labels
    assert len(labels) == 2
    assert 1 in labels
    assert 2 in labels

    # read D
    labels = lb.sweep_label_neighborhood(
        'TATATATATAGCTAGCTAGCTAACTAGCTAGCATCGATCGATCGATC')
    print labels
    assert len(labels) == 1
    assert 3 in labels


def test_label_tag_correctness_save_load():
    lb_pre = LabelHash(20, 1e7, 4)
    filename = utils.get_test_data('test-labels.fa')
    lb_pre.consume_fasta_and_tag_with_labels(filename)

    # save labels to a file
    savepath = utils.get_temp_filename('saved.labels')
    lb_pre.save_labels_and_tags(savepath)

    # trash the old LabelHash
    del lb_pre

    # create new, load labels & tags
    lb = LabelHash(20, 1e7, 4)
    lb.load_labels_and_tags(savepath)

    # read A
    labels = lb.sweep_label_neighborhood(
        'ATCGTGTAAGCTATCGTAATCGTAAGCTCTGCCTAGAGCTAGGCTAGGCTCTGCCTAGAG'
        'CTAGGCTAGGTGTGCTCTGCCTAGAGCTAGGCTAGGTGT')
    print lb.sweep_tag_neighborhood(
        'TTCGTGTAAGCTATCGTAATCGTAAGCTCTGCCTAGAGCTAGGCTAGGCTCTGCCTAGAG'
        'CTAGGCTAGGTGTGCTCTGCTAGAGCTAGGCTAGGTGT')
    print labels
    print len('ATCGTGTAAGCTATCGTAATCGTAAGCTCTGCCTAGAGCTAGGCTAG') - 19
    assert len(labels) == 2
    assert 0 in labels
    assert 1 in labels

    # read B
    labels = lb.sweep_label_neighborhood(
        'GCGTAATCGTAAGCTCTGCCTAGAGCTAGGCTAGCTCTGCCTAGAGCTAGGCTAGGTGTTGGGGATAG'
        'ATAGATAGATGACCTAGAGCTAGGCTAGGTGTTGGGGATAGATAGATAGATGA')
    print labels
    assert len(labels) == 3
    assert 0 in labels
    assert 1 in labels
    assert 2 in labels

    # read C
    labels = lb.sweep_label_neighborhood(
        'TGGGATAGATAGATAGATGACCTAGAGCTAGGCTAGGTGTTGGGGATAGATAGATAGATGACCTAGAG'
        'CTAGGCTAGGTGTTGGGGATAGATAGATAGATGAGTTGGGGATAGATAGATAGATGAGTGTAGATCCA'
        'ACAACACATACA')
    print labels
    assert len(labels) == 2
    assert 1 in labels
    assert 2 in labels

    # read D
    labels = lb.sweep_label_neighborhood(
        'TATATATATAGCTAGCTAGCTAACTAGCTAGCATCGATCGATCGATC')
    print labels
    assert len(labels) == 1
    assert 3 in labels
