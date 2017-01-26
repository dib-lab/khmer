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
from __future__ import print_function, unicode_literals


def print_error(msg):
    """Print the given message to 'stderr'."""
    import sys

    print(msg, file=sys.stderr)


def _split_left_right(name):
    """Split record name at the first whitespace and return both parts.

    RHS is set to an empty string if not present.
    """
    parts = name.split(None, 1)
    lhs, rhs = [parts[0], parts[1] if len(parts) > 1 else '']
    return lhs, rhs


def check_is_pair(record1, record2):
    """Check if the two sequence records belong to the same fragment.

    In an matching pair the records are left and right pairs
    of each other, respectively.  Returns True or False as appropriate.

    Handles both Casava formats: seq/1 and seq/2, and 'seq::... 1::...'
    and 'seq::... 2::...'.

    Also handles the default format of the SRA toolkit's fastq-dump:
    'Accession seq/1'
    """
    if hasattr(record1, 'quality') or hasattr(record2, 'quality'):
        if not (hasattr(record1, 'quality') and hasattr(record2, 'quality')):
            raise ValueError("both records must be same type (FASTA or FASTQ)")

    lhs1, rhs1 = _split_left_right(record1.name)
    lhs2, rhs2 = _split_left_right(record2.name)

    # handle 'name/1'
    if lhs1.endswith('/1') and lhs2.endswith('/2'):
        subpart1 = lhs1.split('/', 1)[0]
        subpart2 = lhs2.split('/', 1)[0]

        if subpart1 and subpart1 == subpart2:
            return True

    # handle '@name 1:rst'
    elif lhs1 == lhs2 and rhs1.startswith('1:') and rhs2.startswith('2:'):
        return True

    # handle @name seq/1
    elif lhs1 == lhs2 and rhs1.endswith('/1') and rhs2.endswith('/2'):
        subpart1 = rhs1.split('/', 1)[0]
        subpart2 = rhs2.split('/', 1)[0]

        if subpart1 and subpart1 == subpart2:
            return True

    return False


def check_is_left(name):
    """Check if the name belongs to a 'left' sequence (/1).

    Returns True or False.

    Handles both Casava formats: seq/1 and 'seq::... 1::...'
    """
    lhs, rhs = _split_left_right(name)
    if lhs.endswith('/1'):              # handle 'name/1'
        return True
    elif rhs.startswith('1:'):          # handle '@name 1:rst'
        return True

    elif rhs.endswith('/1'):            # handles '@name seq/1'
        return True

    return False


def check_is_right(name):
    """Check if the name belongs to a 'right' sequence (/2).

    Returns True or False.

    Handles both Casava formats: seq/2 and 'seq::... 2::...'
    """
    lhs, rhs = _split_left_right(name)
    if lhs.endswith('/2'):              # handle 'name/2'
        return True
    elif rhs.startswith('2:'):          # handle '@name 2:rst'
        return True

    elif rhs.endswith('/2'):            # handles '@name seq/2'
        return True

    return False


class UnpairedReadsError(ValueError):
    """ValueError with refs to the read pair in question."""

    def __init__(self, msg, r1, r2):
        r1_name = "<no read>"
        r2_name = "<no read>"
        if r1:
            r1_name = r1.name
        if r2:
            r2_name = r2.name

        msg = msg + '\n"{0}"\n"{1}"'.format(r1_name, r2_name)
        ValueError.__init__(self, msg)
        self.read1 = r1
        self.read2 = r2


def broken_paired_reader(screed_iter, min_length=None,
                         force_single=False, require_paired=False):
    """Read pairs from a stream.

    A generator that yields singletons and pairs from a stream of FASTA/FASTQ
    records (yielded by 'screed_iter').  Yields (n, is_pair, r1, r2) where
    'r2' is None if is_pair is False.

    The input stream can be fully single-ended reads, interleaved paired-end
    reads, or paired-end reads with orphans, a.k.a. "broken paired".

    Usage::

       for n, is_pair, read1, read2 in broken_paired_reader(...):
          ...

    Note that 'n' behaves like enumerate() and starts at 0, but tracks
    the number of records read from the input stream, so is
    incremented by 2 for a pair of reads.

    If 'min_length' is set, all reads under this length are ignored (even
    if they are pairs).

    If 'force_single' is True, all reads are returned as singletons.
    """
    record = None
    prev_record = None
    num = 0

    if force_single and require_paired:
        raise ValueError("force_single and require_paired cannot both be set!")

    # handle the majority of the stream.
    for record in screed_iter:
        if prev_record:
            if check_is_pair(prev_record, record) and not force_single:
                if min_length and (len(prev_record.sequence) < min_length or
                                   len(record.sequence) < min_length):
                    if require_paired:
                        record = None
                else:
                    yield num, True, prev_record, record  # it's a pair!
                    num += 2
                    record = None
            else:                                   # orphan.
                if require_paired:
                    err = UnpairedReadsError(
                        "Unpaired reads when require_paired is set!",
                        prev_record, record)
                    raise err

                # ignore short reads
                if min_length and len(prev_record.sequence) < min_length:
                    pass
                else:
                    yield num, False, prev_record, None
                    num += 1

        prev_record = record
        record = None

    # handle the last record, if it exists (i.e. last two records not a pair)
    if prev_record:
        if require_paired:
            raise UnpairedReadsError("Unpaired reads when require_paired "
                                     "is set!", prev_record, None)
        if min_length and len(prev_record.sequence) < min_length:
            pass
        else:
            yield num, False, prev_record, None


def write_record(record, fileobj):
    """Write sequence record to 'fileobj' in FASTA/FASTQ format."""
    if hasattr(record, 'quality'):
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

    if hasattr(read1, 'quality'):
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


# vim: set filetype=python tabstop=4 softtabstop=4 shiftwidth=4 expandtab:
# vim: set textwidth=79:
