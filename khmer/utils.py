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
"""Helpful methods for performing common argument-checking tasks in scripts."""
from khmer._oxli.parsing import (check_is_left, check_is_right, check_is_pair,
                                 UnpairedReadsError, _split_left_right,
                                 FastxParser, SplitPairedReader,
                                 BrokenPairedReader)
import itertools


PAIRING_MODES = ('split', 'interleaved', 'single')

def grouper(n, iterable):
    iterable = iter(iterable)
    return iter(lambda: list(itertools.islice(iterable, n)), [])


def print_error(msg):
    """Print the given message to 'stderr'."""
    import sys

    print(msg, file=sys.stderr)


def paired_fastx_handler(samples, pairing_mode, min_length=-1,
                         force_name_match=False, yield_filenames=False, 
                         **kwargs):

    if pairing_mode not in PAIRING_MODES:
        raise ValueError('Pairing mode must be one of {0}'.format(PAIRING_MODES))
    
    if pairing_mode == 'split':
        _samples = grouper(2, samples)
    else:
        _samples = samples

    for group in _samples:
        if pairing_mode == 'split':
            reader = SplitPairedReader(FastxParser(group[0]),
                                       FastxParser(group[1]),
                                       min_length=min_length,
                                       force_name_match=force_name_match)
        elif pairing_mode == 'single':
            reader = BrokenPairedReader(FastxParser(group),
                                        force_single=True,
                                        min_length=min_length,
                                        require_paired=force_name_match)
        else:
            reader = BrokenPairedReader(FastxParser(group),
                                        force_single=False,
                                        min_length=min_length,
                                        require_paired=force_name_match)
        if yield_filenames:
            if pairing_mode == 'split':
                _filename = group[0] + '.pair'
            else:
                _filename = group
            yield _filename, reader
        else:
            yield reader


def write_record(record, fileobj):
    """Write sequence record to 'fileobj' in FASTA/FASTQ format."""
    if hasattr(record, 'quality') and record.quality is not None:
        recstr = '@{name}\n{sequence}\n+\n{quality}\n'.format(
            name=record.name,
            sequence=record.sequence,
            quality=record.quality)
    else:
        recstr = '>{name}\n{sequence}\n'.format(
            name=record.name,
            sequence=record.sequence)

    try:
        fileobj.write(bytes(recstr, 'ascii'))
    except TypeError:
        fileobj.write(recstr)


def write_record_pair(read1, read2, fileobj):
    """Write a pair of sequence records to 'fileobj' in FASTA/FASTQ format."""
    _rec_pair = '@%s\n%s\n+\n%s\n' * 2
    _rec_pair_no_qual = '>%s\n%s\n' * 2

    if hasattr(read1, 'quality') and read1.quality is not None:
        assert hasattr(read2, 'quality')
        recstr = _rec_pair % (read1.name, read1.sequence, read1.quality,
                              read2.name, read2.sequence, read2.quality)

    else:
        recstr = _rec_pair_no_qual % (read1.name, read1.sequence,
                                      read2.name, read2.sequence)

    try:
        fileobj.write(bytes(recstr, 'ascii'))
    except TypeError:
        fileobj.write(recstr)


def clean_input_reads(records):
    """Add a cleaned_seq attribute to records that do not have one

    Use this to convert a stream of records that might not have a
    `cleaned_seq` attribute to one that does. Use this to extend
    Records loaded by `screed.open()`. It is a mistake to apply
    this to a `ReadParser` stream.
    """
    for record in records:
        record.cleaned_seq = record.sequence.upper().replace('N', 'A')
        yield record


class ReadBundle(object):
    def __init__(self, *raw_records):
        self.reads = [i for i in raw_records if i]

    def coverages(self, graph):
        return [graph.get_median_count(r.cleaned_seq)[0] for r in self.reads]

    def coverages_at_least(self, graph, coverage):
        return all(graph.median_at_least(r.cleaned_seq, coverage)
                   for r in self.reads)

    @property
    def num_reads(self):
        return len(self.reads)

    @property
    def total_length(self):
        return sum([len(r.sequence) for r in self.reads])

def grouper(n, iterable):
    iterable = iter(iterable)
    return iter(lambda: list(itertools.islice(iterable, n)), [])


# vim: set filetype=python tabstop=4 softtabstop=4 shiftwidth=4 expandtab:
# vim: set textwidth=79:
