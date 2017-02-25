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
# pylint: disable=missing-docstring,invalid-name

# Tests for the ReadParser and Read classes.
from __future__ import print_function
from __future__ import absolute_import
from khmer import ReadParser, Counttable, Nodegraph
from screed import Record
from . import khmer_tst_utils as utils
import pytest
from functools import reduce  # pylint: disable=redefined-builtin


def test_read_cleaning_consume_seqfile():
    infile = utils.get_test_data('valid-read-testing.fq')

    x = Counttable(15, int(1e6), 4)
    x.consume_seqfile(infile)

    # the relevant read will automatically get uppercased => abundance of 2
    kmer = "caggcgcccaccacc".upper()
    assert x.get(kmer) == 2

    # the 2nd read with this k-mer in it has an N in it.
    kmer = "CCTCATCGGCACCAG"
    assert x.get(kmer) == 2

    # the 2nd read with this k-mer in it has a Z in it
    kmer = "ACTGAGCTTCATGTC"
    assert x.get(kmer) == 2


def test_read_cleaning_consume_read_by_read():
    infile = utils.get_test_data('valid-read-testing.fq')

    x = Counttable(15, int(1e6), 4)
    for read in ReadParser(infile):
        x.consume(read.sequence)          # consume raw sequence

    # the relevant read will be entirely ignored
    # (b/c ReadParser does not uppercase)
    kmer = "caggcgcccaccacc".upper()
    assert x.get(kmer) == 1

    # consume will ignore the invalid base in 2nd read containing this k-mer,
    # so the k-mer will have an abundance of 2.
    kmer = "CCTCATCGGCACCAG"
    assert x.get(kmer) == 2

    # consume will ignore the invalid base in 2nd read containing this k-mer,
    # so the k-mer will have an abundance of 2.
    kmer = "ACTGAGCTTCATGTC"
    assert x.get(kmer) == 2


def test_read_cleaning_consume_read_by_read_cleaned_seq():
    infile = utils.get_test_data('valid-read-testing.fq')

    x = Counttable(15, int(1e6), 4)
    for read in ReadParser(infile):
        x.consume(read.cleaned_seq)       # consume cleaned_seq

    # the relevant read will be cleaned & loaded
    kmer = "caggcgcccaccacc".upper()
    assert x.get(kmer) == 2

    # this k-mer will be correctly loaded
    kmer = "CCTCATCGGCACCAG"
    assert x.get(kmer) == 2

    # this k-mer will be correctly loaded
    kmer = "ACTGAGCTTCATGTC"
    assert x.get(kmer) == 2


def test_read_cleaning_abundance_distribution():
    infile = utils.get_test_data('valid-read-testing.fq')

    x = Counttable(15, int(1e6), 4)
    y = Nodegraph(15, int(1e6), 4)

    x.consume_seqfile(infile)

    dist = x.abundance_distribution(infile, y)
    assert dist[1] == 29                  # k-mers with non-ACGTN => ignored.
    assert dist[2] == 68


def test_read_cleaning_trim_functions_lowercase():
    infile = utils.get_test_data('valid-read-testing.fq')

    # read this in using "approved good" behavior w/cleaned_seq
    x = Counttable(8, int(1e6), 4)
    for read in ReadParser(infile):
        x.consume(read.cleaned_seq)       # consume cleaned_seq

    # all of these functions will fail to do anything, b/c lowercase != valid
    # BUT they will not raise an exception, either.

    s = "caggcgcccaccaccgtgccctccaacctgatggt"
    _, where = x.trim_on_abundance(s, 1)
    assert where == 0

    _, where = x.trim_below_abundance(s, 0)
    print(x.get_kmer_counts(s))
    assert where == 35                    # stays at 35 (abunds all == 0)

    posns = x.find_spectral_error_positions(s, 1)
    assert posns == []


def test_read_cleaning_trim_functions_N():
    infile = utils.get_test_data('valid-read-testing.fq')

    # read this in using "approved good" behavior w/cleaned_seq
    x = Counttable(8, int(1e6), 4)
    for read in ReadParser(infile):
        x.consume(read.cleaned_seq)       # consume cleaned_seq

    s = "ACTGGGCGTAGNCGGTGTCCTCATCGGCACCAGC"
    _, where = x.trim_on_abundance(s, 1)
    assert where == 11

    _, where = x.trim_below_abundance(s, 2)
    assert where == 34

    posns = x.find_spectral_error_positions(s, 1)
    assert posns == [11]


def test_read_cleaning_trim_functions_bad_dna():
    infile = utils.get_test_data('valid-read-testing.fq')

    # read this in using "approved good" behavior w/cleaned_seq
    x = Counttable(8, int(1e6), 4)
    for read in ReadParser(infile):
        x.consume(read.cleaned_seq)       # consume cleaned_seq

    s = "CCGGCGTGGTTZZYAGGTCACTGAGCTTCATGTC"
    _, where = x.trim_on_abundance(s, 1)
    assert where == 11

    _, where = x.trim_below_abundance(s, 2)
    assert where == 34

    posns = x.find_spectral_error_positions(s, 1)
    assert posns == [11]


def test_read_cleaning_output_partitions():
    infile = utils.get_test_data('valid-read-testing.fq')
    savepath = utils.get_temp_filename('foo')

    # read this in using "approved good" behavior w/cleaned_seq
    x = Nodegraph(8, int(1e6), 4)
    for read in ReadParser(infile):
        x.consume(read.cleaned_seq)       # consume cleaned_seq

    kmer = 'caggcgcc'.upper()
    x.add_tag(kmer)
    x.set_partition_id(kmer, 1)

    kmer = 'ACTGGGCG'
    x.add_tag(kmer)
    x.set_partition_id(kmer, 2)

    kmer = 'CCGGCGTG'
    x.add_tag(kmer)
    x.set_partition_id(kmer, 3)

    x.output_partitions(infile, savepath)

    read_names = [read.name for read in ReadParser(savepath)]
    print(read_names)
    assert len(read_names) == 6

    assert '895:1:1:1246:14654 1:N:0:NNNNN\t1' in read_names
    assert '895:1:1:1248:9583 1:N:0:NNNNN\t2' in read_names
    assert '895:1:1:1252:19493 1:N:0:NNNNN\t3' in read_names

    assert 'lowercase_to_uppercase\t1' in read_names
    assert 'n_in_read\t2' in read_names
    assert 'zy_in_read\t3' in read_names


def test_read_cleaning_trim_on_stoptags():
    infile = utils.get_test_data('valid-read-testing.fq')
    savepath = utils.get_temp_filename('foo')

    # read this in using "approved good" behavior w/cleaned_seq
    x = Nodegraph(8, int(1e6), 4)
    for read in ReadParser(infile):
        x.consume(read.cleaned_seq)       # consume cleaned_seq

    # add this as a stop tag
    kmer = 'caggcgcc'.upper()
    x.add_stop_tag(kmer)

    kmer = 'ACTGGGCG'
    x.add_stop_tag(kmer)

    kmer = 'CCGGCGTG'
    x.add_stop_tag(kmer)

    _, pos = x.trim_on_stoptags('caggcgcccaccaccgtgccctccaacctgatggt')
    assert pos == 35                      # no stoptag b/c lowercase => no trim

    _, pos = x.trim_on_stoptags('ACTGGGCGTAGNCGGTGTCCTCATCGGCACCAGC')
    assert pos == 6                       # N ignored

    _, pos = x.trim_on_stoptags('CCGGCGTGGTTZZYAGGTCACTGAGCTTCATGTC')
    assert pos == 6                       # ZZY ignored


# vim: set filetype=python tabstop=4 softtabstop=4 shiftwidth=4 expandtab:
# vim: set textwidth=79:
