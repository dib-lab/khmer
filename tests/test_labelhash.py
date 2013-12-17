#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
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
    expected = [0L, 1L, 2L, 3L]
    for e_label in expected:
        assert e_label in labels
    for a_label in labels:
        assert a_label in expected

def test_get_tag_labels():
    lb = LabelHash(20, 1e7, 4)
    filename = utils.get_test_data('single-read.fq')
    lb.consume_fasta_and_tag_with_labels(filename)
    tag = 173473779682L

    labels = lb.get_tag_labels(tag)
    assert len(labels) == 1
    assert labels.pop() == 0L

def test_consume_partitioned_fasta_and_tag_with_labels():
    lb = LabelHash(20, 1e7, 4)
    filename = utils.get_test_data('real-partition-small.fa')

    total_reads, n_consumed = lb.consume_partitioned_fasta_and_tag_with_labels(filename)
    labels = set()
    for record in screed.open(filename):
        seq = record.sequence
        labels.update(lb.sweep_label_neighborhood(seq, False, False))
    #print lb.n_labels()
    #print labels
    assert len(labels) == 1
    assert labels.pop() == 2L
    assert lb.n_labels() == 1 

def test_consume_fasta_and_tag_with_labels():
    lb = LabelHash(20, 1e7, 4)
    read_1 = 'ACGTAACCGGTTAAACCCGGGTTTAAAACCCCGGGGTTTT'
    filename = utils.get_test_data('test-transcript.fa')

    total_reads, n_consumed = lb.consume_fasta_and_tag_with_labels(filename)

    assert lb.get(read_1[:20])
    assert total_reads == 3
    print lb.n_labels()
    print lb.get_label_dict()
    for tag in lb.get_tagset():
        print tag, khmer.forward_hash(tag, 20)
    for record in screed.open(filename):
        print lb.sweep_tag_neighborhood(record.sequence, 40)
        print lb.sweep_label_neighborhood(record.sequence, 40)
    assert lb.n_labels() == 3

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
    labels = lb.sweep_label_neighborhood('ATCGTGTAAGCTATCGTAATCGTAAGCTCTGCCTAGAGCTAGGCTAGGCTCTGCCTAGAGCTAGGCTAGGTGTGCTCTGCCTAGAGCTAGGCTAGGTGT')
    print lb.sweep_tag_neighborhood('TTCGTGTAAGCTATCGTAATCGTAAGCTCTGCCTAGAGCTAGGCTAGGCTCTGCCTAGAGCTAGGCTAGGTGTGCTCTGCCTAGAGCTAGGCTAGGTGT')
    print labels
    print len('ATCGTGTAAGCTATCGTAATCGTAAGCTCTGCCTAGAGCTAGGCTAG')-19 
    assert len(labels) == 2
    assert 0L in labels
    assert 1L in labels
    
    # read B
    labels = lb.sweep_label_neighborhood('GCGTAATCGTAAGCTCTGCCTAGAGCTAGGCTAGCTCTGCCTAGAGCTAGGCTAGGTGTTGGGGATAGATAGATAGATGACCTAGAGCTAGGCTAGGTGTTGGGGATAGATAGATAGATGA')
    print labels
    assert len(labels) == 3
    assert 0L in labels
    assert 1L in labels
    assert 2L in labels
    
    # read C
    labels = lb.sweep_label_neighborhood('TGGGATAGATAGATAGATGACCTAGAGCTAGGCTAGGTGTTGGGGATAGATAGATAGATGACCTAGAGCTAGGCTAGGTGTTGGGGATAGATAGATAGATGAGTTGGGGATAGATAGATAGATGAGTGTAGATCCAACAACACATACA')
    print labels
    assert len(labels) == 2
    assert 1L in labels
    assert 2L in labels
    
    # read D
    labels = lb.sweep_label_neighborhood('TATATATATAGCTAGCTAGCTAACTAGCTAGCATCGATCGATCGATC')
    print labels
    assert len(labels) == 1
    assert 3L in labels

def test_sweep_tag_neighborhood():
    lb = LabelHash(20, 1e7, 4)
    filename = utils.get_test_data('single-read.fq')
    lb.consume_fasta_and_tag(filename)
    
    tags = lb.sweep_tag_neighborhood('CAGGCGCCCACCACCGTGCCCTCCAACCTGATGGT')
    assert len(tags) == 1
    assert tags.pop() == 173473779682L


def test_sweep_label_neighborhood():
    lb = LabelHash(20, 1e7, 4)
    filename = utils.get_test_data('single-read.fq')
    lb.consume_fasta_and_tag_with_labels(filename)
    
    labels = lb.sweep_label_neighborhood('CAGGCGCCCACCACCGTGCCCTCCAACCTGATGGT')
    assert len(labels) == 1
    assert labels.pop() == 0L
