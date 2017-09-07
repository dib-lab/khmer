# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2011-2015, Michigan State University.
# Copyright (C) 2015, The Regents of the University of California.
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
# pylint: disable=missing-docstring,no-member,protected-access,invalid-name

import khmer

from . import khmer_tst_utils as utils

# Below, 'fakelump.fa' is an artificial data set of 3x1 kb sequences in
# which the last 79 bases are common between the 3 sequences.


def test_fakelump_together():
    fakelump_fa = utils.get_test_data('fakelump.fa')

    ht = khmer.Nodegraph(32, 1e5, 4)
    ht.consume_seqfile_and_tag(fakelump_fa)

    subset = ht.do_subset_partition(0, 0)
    ht.merge_subset(subset)

    (n_partitions, _) = ht.count_partitions()
    assert n_partitions == 1, n_partitions

# try loading stop tags from previously saved


def test_fakelump_stop():
    fakelump_fa = utils.get_test_data('fakelump.fa')
    fakelump_stoptags_txt = utils.get_test_data('fakelump.fa.stoptags.txt')

    ht = khmer.Nodegraph(32, 1e5, 4)
    ht.consume_seqfile_and_tag(fakelump_fa)

    for line in open(fakelump_stoptags_txt):
        ht.add_stop_tag(line.strip())

    subset = ht.do_subset_partition(0, 0, True)
    ht.merge_subset(subset)

    (n_partitions, _) = ht.count_partitions()
    assert n_partitions == 3, n_partitions

# check specific insertion of stop tag


def test_fakelump_stop2():
    fakelump_fa = utils.get_test_data('fakelump.fa')

    ht = khmer.Nodegraph(32, 1e5, 4)
    ht.consume_seqfile_and_tag(fakelump_fa)

    ht.add_stop_tag('GGGGAGGGGTGCAGTTGTGACTTGCTCGAGAG')

    subset = ht.do_subset_partition(0, 0, True)
    ht.merge_subset(subset)

    (n_partitions, _) = ht.count_partitions()
    assert n_partitions == 3, n_partitions

# try repartitioning


def test_fakelump_repartitioning():
    fakelump_fa = utils.get_test_data('fakelump.fa')
    fakelump_fa_foo = utils.get_temp_filename('fakelump.fa.stopfoo')

    ht = khmer.Nodegraph(32, 1e5, 4)
    ht.consume_seqfile_and_tag(fakelump_fa)

    subset = ht.do_subset_partition(0, 0)
    ht.merge_subset(subset)

    (n_partitions, _) = ht.count_partitions()
    assert n_partitions == 1, n_partitions

    # now, break partitions on any k-mer that you see more than once
    # on big excursions, where big excursions are excursions 40 out
    # that encounter more than 82 k-mers.  This should specifically
    # identify our connected sequences in fakelump...

    EXCURSION_DISTANCE = 40
    EXCURSION_KMER_THRESHOLD = 82
    EXCURSION_KMER_COUNT_THRESHOLD = 1
    counting = khmer.Countgraph(32, 1e5, 4)

    ht.repartition_largest_partition(counting,
                                     EXCURSION_DISTANCE,
                                     EXCURSION_KMER_THRESHOLD,
                                     EXCURSION_KMER_COUNT_THRESHOLD)

    ht.save_stop_tags(fakelump_fa_foo)

    # ok, now re-do everything with these stop tags, specifically.

    ht = khmer.Nodegraph(32, 1e5, 4)
    ht.consume_seqfile_and_tag(fakelump_fa)
    ht.load_stop_tags(fakelump_fa_foo)

    subset = ht.do_subset_partition(0, 0, True)
    ht.merge_subset(subset)

    (n_partitions, _) = ht.count_partitions()
    assert n_partitions == 6, n_partitions


def test_fakelump_load_stop_tags_trunc():
    fakelump_fa = utils.get_test_data('fakelump.fa')
    fakelump_fa_foo = utils.get_temp_filename('fakelump.fa.stopfoo')

    ht = khmer.Nodegraph(32, 1e5, 4)
    ht.consume_seqfile_and_tag(fakelump_fa)

    subset = ht.do_subset_partition(0, 0)
    ht.merge_subset(subset)

    (n_partitions, _) = ht.count_partitions()
    assert n_partitions == 1, n_partitions

    # now, break partitions on any k-mer that you see more than once
    # on big excursions, where big excursions are excursions 40 out
    # that encounter more than 82 k-mers.  This should specifically
    # identify our connected sequences in fakelump...

    EXCURSION_DISTANCE = 40
    EXCURSION_KMER_THRESHOLD = 82
    EXCURSION_KMER_COUNT_THRESHOLD = 1
    counting = khmer.Countgraph(32, 1, 1, primes=[5, 7, 11, 13])

    ht.repartition_largest_partition(counting,
                                     EXCURSION_DISTANCE,
                                     EXCURSION_KMER_THRESHOLD,
                                     EXCURSION_KMER_COUNT_THRESHOLD)

    ht.save_stop_tags(fakelump_fa_foo)
    data = open(fakelump_fa_foo, 'rb').read()

    fp = open(fakelump_fa_foo, 'wb')
    fp.write(data[:10])
    fp.close()

    # ok, now try loading these stop tags; should fail.
    ht = khmer.Nodegraph(32, 1, 1, primes=[5, 7, 11, 13])
    ht.consume_seqfile_and_tag(fakelump_fa)

    try:
        ht.load_stop_tags(fakelump_fa_foo)
        assert 0, "this test should fail"
    except OSError:
        pass


def test_fakelump_load_stop_tags_notexist():
    fakelump_fa_foo = utils.get_temp_filename('fakelump.fa.stopfoo')

    # ok, now try loading these stop tags; should fail.
    ht = khmer.Nodegraph(32, 1, 1, primes=[5, 7, 11, 13])

    try:
        ht.load_stop_tags(fakelump_fa_foo)
        assert 0, "this test should fail"
    except OSError:
        pass
