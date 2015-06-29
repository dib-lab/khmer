from __future__ import print_function, unicode_literals
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# Convenience functions for performing common argument-checking tasks in
# scripts.


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

        assert subpart1
        if subpart1 == subpart2:
            return True

    # handle '@name 1:rst'
    elif lhs1 == lhs2 and rhs1.startswith('1:') and rhs2.startswith('2:'):
        return True

    # handle @name seq/1
    elif lhs1 == lhs2 and rhs1.endswith('/1') and rhs2.endswith('/2'):
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
    n = 0

    if force_single and require_paired:
        raise ValueError("force_single and require_paired cannot both be set!")

    # handle the majority of the stream.
    for record in screed_iter:
        # ignore short reads
        if min_length and len(record.sequence) < min_length:
            record = None
            continue

        if prev_record:
            if check_is_pair(prev_record, record) and not force_single:
                yield n, True, prev_record, record  # it's a pair!
                n += 2
                record = None
            else:                                   # orphan.
                if require_paired:
                    raise ValueError("Unpaired reads when require_paired"
                                     " is set!")
                yield n, False, prev_record, None
                n += 1

        prev_record = record
        record = None

    # handle the last record, if it exists (i.e. last two records not a pair)
    if prev_record:
        if require_paired:
            raise ValueError("Unpaired reads when require_paired is set!")
        yield n, False, prev_record, None


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
        fileobj.write(bytes(recstr, 'utf-8'))
    except TypeError:
        fileobj.write(recstr)


def write_record_pair(read1, read2, fileobj):
    """Write a pair of sequence records to 'fileobj' in FASTA/FASTQ format."""
    if hasattr(read1, 'quality'):
        assert hasattr(read2, 'quality')
    write_record(read1, fileobj)
    write_record(read2, fileobj)


# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
