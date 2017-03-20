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
from khmer import _Countgraph, _Counttable, _SmallCountgraph, _SmallCounttable
from khmer import _Nodegraph, _Nodetable
from khmer import ReadParser, Nodegraph, _GraphLabels
from . import khmer_tst_utils as utils
import pytest


PRIMES_1m = [1000003, 1009837]


# all the table types!
@pytest.fixture(params=[_Countgraph, _Counttable, _SmallCountgraph,
                        _SmallCounttable, _Nodegraph, _Nodetable])
def Tabletype(request):
    return request.param


# all the counting types!
@pytest.fixture(params=[_Countgraph, _Counttable, _SmallCountgraph,
                        _SmallCounttable])
def Countingtype(request):
    return request.param


# all the graph types!
@pytest.fixture(params=[_Countgraph, _Nodegraph])
def Graphtype(request):
    return request.param


@pytest.yield_fixture
def reads():
    infile = utils.get_test_data('valid-read-testing.fq')
    reads = ReadParser(infile)
    yield reads
    reads.close()


def test_read_cleaning_consume_seqfile(Countingtype):
    infile = utils.get_test_data('valid-read-testing.fq')

    x = Countingtype(15, PRIMES_1m)
    x.consume_seqfile(infile)

    # the relevant read will automatically get uppercased => abundance of 2
    kmer = "caggcgcccaccacc".upper()
    assert x.get(kmer) == 2

    # the 2nd read with this k-mer in it has an N in it.
    kmer = "CCTCATCGGCACCAG"
    assert x.get(kmer) == 1               # this should be 2 in the future

    # the 2nd read with this k-mer in it has a Z in it
    kmer = "ACTGAGCTTCATGTC"
    assert x.get(kmer) == 1               # this should be 2 in the future


def test_read_cleaning_consume_read_by_read(Countingtype, reads):
    x = Countingtype(15, PRIMES_1m)
    for read in reads:
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


def test_read_cleaning_consume_read_by_read_cleaned_seq(Countingtype, reads):
    x = Countingtype(15, PRIMES_1m)
    for read in reads:
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


def test_read_cleaning_abundance_distribution(Countingtype):
    infile = utils.get_test_data('valid-read-testing.fq')

    x = Countingtype(15, PRIMES_1m)
    y = _Nodegraph(15, PRIMES_1m)

    x.consume_seqfile(infile)

    dist = x.abundance_distribution(infile, y)
    assert dist[1] == 41
    assert dist[2] == 42


def test_read_cleaning_trim_functions_lowercase(Tabletype, reads):
    # read this in using "approved good" behavior w/cleaned_seq
    x = Tabletype(8, PRIMES_1m)
    for read in reads:
        x.consume(read.cleaned_seq)       # consume cleaned_seq

    s = "caggcgcccaccaccgtgccctccaacctgatggt"
    _, where = x.trim_on_abundance(s, 1)
    assert where == 35                    # in future, should be '0'

    _, where = x.trim_below_abundance(s, 2)
    assert where == 35                    # in future, should be ?? @CTB

    posns = x.find_spectral_error_positions(s, 1)
    assert posns == []                    # in future, should do same


def test_read_cleaning_trim_functions_N(Tabletype, reads):
    # read this in using "approved good" behavior w/cleaned_seq
    x = Tabletype(8, PRIMES_1m)
    for read in reads:
        x.consume(read.cleaned_seq)       # consume cleaned_seq

    s = "ACTGGGCGTAGNCGGTGTCCTCATCGGCACCAGC"
    _, where = x.trim_on_abundance(s, 1)
    assert where == 0                     # in future, should be ??

    _, where = x.trim_below_abundance(s, 2)
    assert where == 0                     # in future, should be ??

    with pytest.raises(ValueError):
        posns = x.find_spectral_error_positions(s, 1)
    # assert posns == []                    # in future, should return []


def test_read_cleaning_trim_functions_bad_dna(Tabletype, reads):
    # read this in using "approved good" behavior w/cleaned_seq
    x = Tabletype(8, PRIMES_1m)
    for read in reads:
        x.consume(read.cleaned_seq)       # consume cleaned_seq

    s = "CCGGCGTGGTTZZYAGGTCACTGAGCTTCATGTC"
    _, where = x.trim_on_abundance(s, 1)
    assert where == 0                     # in future, should be ??

    _, where = x.trim_below_abundance(s, 2)
    assert where == 0                     # in future, should be ??

    with pytest.raises(ValueError):
        posns = x.find_spectral_error_positions(s, 1)
    # assert posns == []                    # in future, should return [11]


def test_read_cleaning_output_partitions(Graphtype):
    infile = utils.get_test_data('valid-read-testing.fq')
    savepath = utils.get_temp_filename('foo')

    # read this in using "approved good" behavior w/cleaned_seq
    x = Graphtype(8, PRIMES_1m)
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
    assert len(read_names) == 4

    print(read_names)
    assert '895:1:1:1246:14654 1:N:0:NNNNN\t1\t1' in read_names
    assert '895:1:1:1248:9583 1:N:0:NNNNN\t2\t2' in read_names
    assert '895:1:1:1252:19493 1:N:0:NNNNN\t3\t3' in read_names

    assert 'lowercase_to_uppercase\t5\t1' in read_names

    assert 'n_in_read\t6\t2' not in read_names
    assert 'zy_in_read\t7\t3' not in read_names


def test_read_cleaning_trim_on_stoptags(Graphtype):
    infile = utils.get_test_data('valid-read-testing.fq')

    # read this in using "approved good" behavior w/cleaned_seq
    x = Graphtype(8, PRIMES_1m)
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
    assert pos == 6                       # should be ?? in future

    _, pos = x.trim_on_stoptags('CCGGCGTGGTTZZYAGGTCACTGAGCTTCATGTC')
    assert pos == 0                       # should be 6 in future

    _, pos = x.trim_on_stoptags('CCGGCGTGGTTZZYAGGTCACTGAGCTTCATGTC')
    assert pos == 0                       # should be 6 in future


def test_consume_seqfile_and_tag(Graphtype):
    infile = utils.get_test_data('valid-read-testing.fq')

    # read this in consume_and_tag
    x = Graphtype(8, PRIMES_1m)
    x.consume_seqfile_and_tag(infile)
    _, n_tags = x.count_partitions()
    assert n_tags == 4                    # total # of tags


def test_consume_partitioned_seqfile(Graphtype):
    infile = utils.get_test_data('valid-read-testing.fq')

    # read this in consume_and_tag
    x = Graphtype(15, PRIMES_1m)
    x.consume_partitioned_fasta(infile)
    n_partitions, n_tags = x.count_partitions()
    assert n_partitions == 4
    assert n_tags == 0


def test_output_partitioned_file(Graphtype):
    infile = utils.get_test_data('valid-read-testing.fq')
    savepath = utils.get_temp_filename('foo')

    # read this in consume_and_tag
    x = Graphtype(15, PRIMES_1m)
    x.consume_partitioned_fasta(infile)
    x.output_partitions(infile, savepath)

    read_names = [read.name for read in ReadParser(savepath)]
    read_names = set(read_names)

    good_names = ['895:1:1:1246:14654 1:N:0:NNNNN\t1\t5',
                  '895:1:1:1248:9583 1:N:0:NNNNN\t2\t2',
                  '895:1:1:1252:19493 1:N:0:NNNNN\t3\t3',
                  '895:1:1:1255:18861 1:N:0:NNNNN\t4\t8',
                  'lowercase_to_uppercase\t5\t5',
                  '895:1:1:1255:18861 1:N:0:NNNNN\t8\t8']
    good_names = set(good_names)

    assert good_names == read_names


def test_consume_seqfile_and_tag_with_labels(Graphtype):
    infile = utils.get_test_data('valid-read-testing.fq')

    # read this in consume_and_tag
    graph = Graphtype(15, PRIMES_1m)
    x = _GraphLabels(graph)
    x.consume_seqfile_and_tag_with_labels(infile)

    assert x.n_labels() == 6


def test_consume_partitioned_seqfile_and_label(Graphtype):
    infile = utils.get_test_data('valid-read-testing.fq')

    # read this in consume_and_tag
    graph = Graphtype(15, PRIMES_1m)
    x = _GraphLabels(graph)
    x.consume_partitioned_fasta_and_tag_with_labels(infile)

    assert x.n_labels() == 6


# vim: set filetype=python tabstop=4 softtabstop=4 shiftwidth=4 expandtab:
# vim: set textwidth=79:
