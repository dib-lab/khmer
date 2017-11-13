# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2013-2015, Michigan State University.
# Copyright (C) 2015-2016, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org
# pylint: disable=missing-docstring,protected-access,no-member,invalid-name

import os
import khmer
from khmer import GraphLabels, Nodegraph, Countgraph
import screed

import pytest

from . import khmer_tst_utils as utils


def teardown():
    utils.cleanup()

#
# @camillescott TODO: more tests!
#  * thread-safety


@pytest.mark.huge
def test_toobig():
    try:
        GraphLabels.NodeGraphLabels(20, 1e13, 1)
        assert 0, "This should fail."
    except MemoryError as err:
        print(str(err))


def test_error_create():
    try:
        GraphLabels.NodeGraphLabels(None, None, None)
        assert 0, "This should fail."
    except TypeError as err:
        print(str(err))


def test_n_labels():
    lh = GraphLabels.NodeGraphLabels(20, 1e7, 4)
    filename = utils.get_test_data('test-labels.fa')
    lh.consume_seqfile_and_tag_with_labels(filename)

    print(lh.n_labels)
    assert lh.n_labels == 4


def test_get_all_labels():
    lb = GraphLabels.NodeGraphLabels(20, 1e7, 4)
    filename = utils.get_test_data('test-labels.fa')
    lb.consume_seqfile_and_tag_with_labels(filename)

    labels = list(lb.labels())
    expected = [0, 1, 2, 3]
    for e_label in expected:
        assert e_label in labels
    for a_label in labels:
        assert a_label in expected


def test_get_labels_save_load():
    lb_pre = GraphLabels.NodeGraphLabels(20, 1e7, 4)
    filename = utils.get_test_data('test-labels.fa')
    lb_pre.consume_seqfile_and_tag_with_labels(filename)

    # save labels to a file
    savepath = utils.get_temp_filename('saved.labels')
    lb_pre.save_labels_and_tags(savepath)

    # trash the old GraphLabels
    del lb_pre

    # create new, load labels & tags
    graph = Nodegraph(20, 1e7, 4)
    lb = GraphLabels.load(savepath, graph)

    labels = list(lb.labels())
    expected = [0, 1, 2, 3]
    for e_label in expected:
        assert e_label in labels
    for a_label in labels:
        assert a_label in expected


def test_get_labels_save_load_wrong_ksize():
    lb_pre = GraphLabels.NodeGraphLabels(19, 1e7, 4)
    filename = utils.get_test_data('test-labels.fa')
    lb_pre.consume_seqfile_and_tag_with_labels(filename)

    # save labels to a file
    savepath = utils.get_temp_filename('saved.labels')
    lb_pre.save_labels_and_tags(savepath)

    # trash the old GraphLabels
    del lb_pre

    # create new, load labels & tags
    lb = GraphLabels.NodeGraphLabels(20, 1e7, 4)
    try:
        lb.load_labels_and_tags(savepath)
        assert 0, "this should not succeed - different ksize"
    except OSError as err:
        print(str(err))
        assert "Incorrect k-mer size 19" in str(err)


def test_save_load_corrupted():
    lb_pre = GraphLabels.NodeGraphLabels(20, 1e7, 4)
    filename = utils.get_test_data('test-labels.fa')
    lb_pre.consume_seqfile_and_tag_with_labels(filename)

    # save labels to a file
    savepath = utils.get_temp_filename('saved.labels')
    lb_pre.save_labels_and_tags(savepath)

    # trash the old GraphLabels
    del lb_pre

    lb = GraphLabels.NodeGraphLabels(20, 1e7, 4)

    # produce all possible truncated versions of this file
    data = open(savepath, 'rb').read()
    for i in range(len(data)):
        truncated = utils.get_temp_filename('trunc.labels')
        fp = open(truncated, 'wb')
        fp.write(data[:i])
        fp.close()

        try:
            lb.load_labels_and_tags(truncated)
            assert 0, "this should not succeed -- truncated file len %d" % (i,)
        except OSError as err:
            print('expected failure for', i, ': ', str(err))


# note: if run as root, will fail b/c root can write to anything
@pytest.mark.noroot
def test_save_fail_readonly():
    lb_pre = GraphLabels.NodeGraphLabels(20, 1e7, 4)
    filename = utils.get_test_data('test-labels.fa')
    lb_pre.consume_seqfile_and_tag_with_labels(filename)

    # save labels to a file
    savepath = utils.get_temp_filename('saved.labels')
    fp = open(savepath, 'w')
    fp.close()

    os.chmod(savepath, 0x444)

    try:
        lb_pre.save_labels_and_tags(savepath)
        assert 0, "this should fail: read-only file"
    except OSError as err:
        print(str(err))


def test_get_tag_labels():
    lb = GraphLabels.NodeGraphLabels(20, 1e7, 4)
    filename = utils.get_test_data('single-read.fq')
    lb.consume_seqfile_and_tag_with_labels(filename)
    tag = 173473779682

    labels = list(lb.get_tag_labels(tag))
    assert len(labels) == 1
    assert labels.pop() == 0


def test_get_labels_for_sequence():
    lb = GraphLabels.NodeGraphLabels(20, 1e7, 4)
    filename = utils.get_test_data('single-read.fq')
    lb.consume_seqfile_and_tag_with_labels(filename)

    seq = [r.sequence for r in screed.open(filename)][0]
    labels = list(lb.get_labels_for_sequence(seq))

    tag = 173473779682
    labels2 = list(lb.get_tag_labels(tag))

    assert labels == labels2
    assert len(labels) == 1
    assert labels.pop() == 0


def test_link_tag_and_label():
    lb = GraphLabels.NodeGraphLabels(20, 1, 1)

    tag = 173473779682
    lb.add_tag(tag)
    lb.link_tag_and_label(tag, 1)

    labels = list(lb.get_tag_labels(tag))
    assert len(labels) == 1
    assert labels.pop() == 1


def test_link_tag_and_label_using_string():
    lb = GraphLabels.NodeGraphLabels(20, 1, 1)

    kmer = lb.graph.reverse_hash(173473779682)
    lb.add_tag(kmer)
    lb.link_tag_and_label(kmer, 1)

    labels = list(lb.get_tag_labels(kmer))
    assert len(labels) == 1
    assert labels.pop() == 1


def test_link_tag_and_label_using_string_2():
    lb = GraphLabels.NodeGraphLabels(20, 1, 1)

    tag = 173473779682
    kmer = lb.graph.reverse_hash(tag)
    lb.add_tag(kmer)
    lb.link_tag_and_label(kmer, 1)

    labels = list(lb.get_tag_labels(tag))    # <-- use 'tag' instead of 'kmer'
    assert len(labels) == 1
    assert labels.pop() == 1


def test_consume_seqfile_and_tag_with_labels():
    lb = GraphLabels.NodeGraphLabels(20, 1e7, 4)
    read_1 = 'ACGTAACCGGTTAAACCCGGGTTTAAAACCCCGGGGTTTT'
    filename = utils.get_test_data('test-transcript.fa')

    total_reads, _ = lb.consume_seqfile_and_tag_with_labels(filename)
    print("doing get")
    assert lb.graph.get(read_1[:20])
    assert total_reads == 3
    print("doing n_labels")
    print(lb.n_labels)
    print("doing all labels")
    print(lb.labels())
    print("get tagset")
    for tag in lb.graph.get_tagset():
        print("forward hash")
        print(tag, khmer.forward_hash(tag, 20))
    for record in screed.open(filename):
        print("Sweeping tags")
        print(lb.sweep_tag_neighborhood(record.sequence, 40))
        print("Sweeping labels...")
        print(lb.sweep_label_neighborhood(record.sequence, 40))
    assert lb.n_labels == 3


def test_consume_partitioned_fasta_and_tag_with_labels():
    lb = GraphLabels.NodeGraphLabels(20, 1e7, 4)
    filename = utils.get_test_data('real-partition-small.fa')

    lb.consume_partitioned_fasta_and_tag_with_labels(
        filename)
    labels = set()
    for record in screed.open(filename):
        seq = record.sequence
        labels.update(lb.sweep_label_neighborhood(seq, 0, False, False))
    # print(lb.n_labels())
    # print(labels)
    assert len(labels) == 1
    assert labels.pop() == 2
    assert lb.n_labels == 1


def test_consume_sequence_and_tag_with_labels():
    lb = GraphLabels.NodeGraphLabels(20, 1e6, 4)
    label = 0
    sequence = 'ATGCATCGATCGATCGATCGATCGATCGATCGATCGATCG'

    lb.consume_sequence_and_tag_with_labels(sequence, label)
    labels = set()
    labels.update(lb.sweep_label_neighborhood(sequence))

    assert label in labels
    assert len(labels) == 1


def test_consume_sequence_and_tag_with_labels_2():
    lb = GraphLabels.NodeGraphLabels(20, 1e6, 4)
    label = 56                            # randomly chosen / non-zero
    sequence = 'ATGCATCGATCGATCGATCGATCGATCGATCGATCGATCG'

    lb.consume_sequence_and_tag_with_labels(sequence, label)
    labels = set()
    labels.update(lb.sweep_label_neighborhood(sequence))

    assert label in labels
    assert len(labels) == 1


def test_sweep_tag_neighborhood():
    lb = GraphLabels.NodeGraphLabels(20, 1e7, 4)
    filename = utils.get_test_data('single-read.fq')
    lb.graph.consume_seqfile_and_tag(filename)

    tags = lb.sweep_tag_neighborhood('CAGGCGCCCACCACCGTGCCCTCCAACCTGATGGT')
    assert len(tags) == 1
    assert list(tags) == [173473779682]


def test_sweep_label_neighborhood():
    lb = GraphLabels.NodeGraphLabels(20, 1e7, 4)
    filename = utils.get_test_data('single-read.fq')
    lb.consume_seqfile_and_tag_with_labels(filename)

    labels = list(lb.sweep_label_neighborhood(
        'CAGGCGCCCACCACCGTGCCCTCCAACCTGATGGT'))
    assert len(labels) == 1
    assert labels.pop() == 0

#
# * The test data set as four reads: A, B, C, and D
# * Overlaps are A <-> B <-> C, with D on its own
# * Thus, traversing from A should find labels from A and B,
#   traversing from B should find labels from A, B, and C,
#   and traversing from C should find labels from B and C


def test_label_tag_correctness():
    lb = GraphLabels.NodeGraphLabels(20, 1e7, 4)
    filename = utils.get_test_data('test-labels.fa')
    lb.consume_seqfile_and_tag_with_labels(filename)

    # read A
    labels = list(lb.sweep_label_neighborhood(
        'ATCGTGTAAGCTATCGTAATCGTAAGCTCTGCCTAGAGCTAGGCTAGGCTCTGCCTAGAG'
        'CTAGGCTAGGTGTGCTCTGCCTAGAGCTAGGCTAGGTGT'))
    print(lb.sweep_tag_neighborhood(
        'TTCGTGTAAGCTATCGTAATCGTAAGCTCTGCCTAGAGCTAGGCTAGGCTCTGCCTAGAG'
        'CTAGGCTAGGTGTGCTCTGCTAGAGCTAGGCTAGGTGT'))
    print(labels)
    print(len('ATCGTGTAAGCTATCGTAATCGTAAGCTCTGCCTAGAGCTAGGCTAG') - 19)
    assert len(labels) == 2
    assert 0 in labels
    assert 1 in labels

    # read B
    labels = list(lb.sweep_label_neighborhood(
        'GCGTAATCGTAAGCTCTGCCTAGAGCTAGGCTAGCTCTGCCTAGAGCTAGGCTAGGTGTTGGGGATAG'
        'ATAGATAGATGACCTAGAGCTAGGCTAGGTGTTGGGGATAGATAGATAGATGA'))
    print(labels)
    assert len(labels) == 3
    assert 0 in labels
    assert 1 in labels
    assert 2 in labels

    # read C
    labels = list(lb.sweep_label_neighborhood(
        'TGGGATAGATAGATAGATGACCTAGAGCTAGGCTAGGTGTTGGGGATAGATAGATAGATGACCTAGAG'
        'CTAGGCTAGGTGTTGGGGATAGATAGATAGATGAGTTGGGGATAGATAGATAGATGAGTGTAGATCCA'
        'ACAACACATACA'))
    print(labels)
    assert len(labels) == 2
    assert 1 in labels
    assert 2 in labels

    # read D
    labels = list(lb.sweep_label_neighborhood(
        'TATATATATAGCTAGCTAGCTAACTAGCTAGCATCGATCGATCGATC'))
    print(labels)
    assert len(labels) == 1
    assert 3 in labels


def test_counting_label_tag_correctness():
    lb = GraphLabels.CountGraphLabels(20, 1e7, 4)
    filename = utils.get_test_data('test-labels.fa')
    lb.consume_seqfile_and_tag_with_labels(filename)

    # read A
    labels = list(lb.sweep_label_neighborhood(
        'ATCGTGTAAGCTATCGTAATCGTAAGCTCTGCCTAGAGCTAGGCTAGGCTCTGCCTAGAG'
        'CTAGGCTAGGTGTGCTCTGCCTAGAGCTAGGCTAGGTGT'))
    print(lb.sweep_tag_neighborhood(
        'TTCGTGTAAGCTATCGTAATCGTAAGCTCTGCCTAGAGCTAGGCTAGGCTCTGCCTAGAG'
        'CTAGGCTAGGTGTGCTCTGCTAGAGCTAGGCTAGGTGT'))
    print(labels)
    print(len('ATCGTGTAAGCTATCGTAATCGTAAGCTCTGCCTAGAGCTAGGCTAG') - 19)
    assert len(labels) == 2
    assert 0 in labels
    assert 1 in labels

    # read B
    labels = list(lb.sweep_label_neighborhood(
        'GCGTAATCGTAAGCTCTGCCTAGAGCTAGGCTAGCTCTGCCTAGAGCTAGGCTAGGTGTTGGGGATAG'
        'ATAGATAGATGACCTAGAGCTAGGCTAGGTGTTGGGGATAGATAGATAGATGA'))
    print(labels)
    assert len(labels) == 3
    assert 0 in labels
    assert 1 in labels
    assert 2 in labels

    # read C
    labels = list(lb.sweep_label_neighborhood(
        'TGGGATAGATAGATAGATGACCTAGAGCTAGGCTAGGTGTTGGGGATAGATAGATAGATGACCTAGAG'
        'CTAGGCTAGGTGTTGGGGATAGATAGATAGATGAGTTGGGGATAGATAGATAGATGAGTGTAGATCCA'
        'ACAACACATACA'))
    print(labels)
    assert len(labels) == 2
    assert 1 in labels
    assert 2 in labels

    # read D
    labels = list(lb.sweep_label_neighborhood(
        'TATATATATAGCTAGCTAGCTAACTAGCTAGCATCGATCGATCGATC'))
    print(labels)
    assert len(labels) == 1
    assert 3 in labels


def test_label_tag_correctness_save_load():
    lb_pre = GraphLabels.NodeGraphLabels(20, 1e7, 4)
    filename = utils.get_test_data('test-labels.fa')
    lb_pre.consume_seqfile_and_tag_with_labels(filename)

    # save labels to a file
    savepath = utils.get_temp_filename('saved.labels')
    lb_pre.save_labels_and_tags(savepath)

    # trash the old GraphLabels
    del lb_pre

    # create new, load labels & tags
    lb = GraphLabels.NodeGraphLabels(20, 1e7, 4)
    lb.load_labels_and_tags(savepath)

    # read A
    labels = list(lb.sweep_label_neighborhood(
        'ATCGTGTAAGCTATCGTAATCGTAAGCTCTGCCTAGAGCTAGGCTAGGCTCTGCCTAGAG'
        'CTAGGCTAGGTGTGCTCTGCCTAGAGCTAGGCTAGGTGT'))
    print(lb.sweep_tag_neighborhood(
        'TTCGTGTAAGCTATCGTAATCGTAAGCTCTGCCTAGAGCTAGGCTAGGCTCTGCCTAGAG'
        'CTAGGCTAGGTGTGCTCTGCTAGAGCTAGGCTAGGTGT'))
    print(labels)
    print(len('ATCGTGTAAGCTATCGTAATCGTAAGCTCTGCCTAGAGCTAGGCTAG') - 19)
    assert len(labels) == 2
    assert 0 in labels
    assert 1 in labels

    # read B
    labels = list(lb.sweep_label_neighborhood(
        'GCGTAATCGTAAGCTCTGCCTAGAGCTAGGCTAGCTCTGCCTAGAGCTAGGCTAGGTGTTGGGGATAG'
        'ATAGATAGATGACCTAGAGCTAGGCTAGGTGTTGGGGATAGATAGATAGATGA'))
    print(labels)
    assert len(labels) == 3
    assert 0 in labels
    assert 1 in labels
    assert 2 in labels

    # read C
    labels = list(lb.sweep_label_neighborhood(
        'TGGGATAGATAGATAGATGACCTAGAGCTAGGCTAGGTGTTGGGGATAGATAGATAGATGACCTAGAG'
        'CTAGGCTAGGTGTTGGGGATAGATAGATAGATGAGTTGGGGATAGATAGATAGATGAGTGTAGATCCA'
        'ACAACACATACA'))
    print(labels)
    assert len(labels) == 2
    assert 1 in labels
    assert 2 in labels

    # read D
    labels = list(lb.sweep_label_neighborhood(
        'TATATATATAGCTAGCTAGCTAACTAGCTAGCATCGATCGATCGATC'))
    print(labels)
    assert len(labels) == 1
    assert 3 in labels


def test_load_wrong_filetype():
    lb = GraphLabels.NodeGraphLabels(20, 1e7, 4)

    # try to load a tagset
    filename = utils.get_test_data('goodversion-k32.tagset')
    try:
        lb.load_labels_and_tags(filename)
        assert 0, "this should not succeed - bad file type"
    except OSError as err:
        print(str(err))
        assert "Incorrect file format type" in str(err)

    # try to load a nonsense file
    filename = utils.get_test_data('all-A.fa')
    try:
        lb.load_labels_and_tags(filename)
        assert 0, "this should not succeed - bad file signature"
    except OSError as err:
        print(str(err))
        assert "Incorrect file signature" in str(err)


def test_load_wrong_fileversion():
    lb = GraphLabels.NodeGraphLabels(20, 1e7, 4)

    # try to load a tagset from an old version
    filename = utils.get_test_data('badversion-k32.tagset')
    try:
        lb.load_labels_and_tags(filename)
        assert 0, "this should not succeed - bad file type"
    except OSError as err:
        print(str(err))
        assert "Incorrect file format version" in str(err)
