#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# Convenience functions for performing common argument-checking tasks in
# scripts.


def print_error(msg):
    """
    Prints the given message to 'stderr'.
    """

    import sys

    print >>sys.stderr, msg


def check_is_pair(record1, record2):
    """
    Checks to see if the two sequence records are left and right pairs
    of each other, respectively.  Returns True or False as appropriate.

    Handles both Casava formats: seq/1 and seq/2, and 'seq::... 1::...'
    and 'seq::... 2::...'.
    """
    if hasattr(record1, 'accuracy') or hasattr(record2, 'accuracy'):
        if not (hasattr(record1, 'accuracy') and hasattr(record2, 'accuracy')):
            raise ValueError("both records must be same type (FASTA or FASTQ)")

    name1 = record1.name
    name2 = record2.name

    if ' ' in name1:                        # handle '@name 1:rst'
        name1, rest1 = record1.name.split(' ', 1)
        name2, rest2 = record2.name.split(' ', 1)

        if name1 == name2 and \
           rest1.startswith('1:') and rest2.startswith('2:'):
            return True

    elif name1.endswith('/1') and name2.endswith('/2'):  # handle name/1
        subpart1 = name1.split('/', 1)[0]
        subpart2 = name2.split('/', 1)[0]

        assert subpart1
        if subpart1 == subpart2:
            return True

    return False


def broken_paired_reader(screed_iter, min_length=None, force_single=False):
    """
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
                yield n, False, prev_record, None
                n += 1

        prev_record = record
        record = None

    # handle the last record, if it exists (i.e. last two records not a pair)
    if prev_record:
        yield n, False, prev_record, None


def write_record(record, fileobj):
    """
    Writes sequence record to 'fileobj' in FASTA/FASTQ format.
    """
    if hasattr(record, 'accuracy'):
        fileobj.write(
            '@{name}\n{seq}\n'
            '+\n{acc}\n'.format(name=record.name,
                                seq=record.sequence,
                                acc=record.accuracy))
    else:
        fileobj.write(
            '>{name}\n{seq}\n'.format(name=record.name,
                                      seq=record.sequence))


def write_record_pair(read1, read2, fileobj):
    """
    Writes pair of sequence records to 'fileobj' in FASTA/FASTQ format.
    """
    if hasattr(read1, 'accuracy'):
        assert hasattr(read2, 'accuracy')
    write_record(read1, fileobj)
    write_record(read2, fileobj)


# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
