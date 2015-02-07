#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
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
    is_fastq = False
    if hasattr(record1, 'accuracy'):
        is_fastq = True
        assert hasattr(record2, 'accuracy')

    name1 = record1.name
    name2 = record2.name

    if is_fastq and ' ' in name1:                        # handle '@name 1:rst'
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


def broken_paired_reader(screed_iter):
    record = None
    prev_record = None

    for n, record in enumerate(screed_iter):
        if prev_record:
            if check_is_pair(prev_record, record):
                yield n, True, prev_record, record  # it's a pair!
                record = None
            else:                                   # orphan.
                yield n, False, prev_record, None

        prev_record = record

    if prev_record:
        if check_is_pair(prev_record, record):
            yield n, True, prev_record, record  # it's a pair!
            record = None
        else:                                   # orphan.
            yield n, False, prev_record, None

    if record:                     # guaranteed to be orphan
        n += 1
        yield n, False, record, None


def write_record(record, fp):
    """
    Writes sequence record to 'fp' in FASTA/FASTQ format.
    """
    if hasattr(record, 'accuracy'):
        fp.write(
            '@{name}\n{seq}\n'
            '+\n{acc}\n'.format(name=record.name,
                                seq=record.sequence,
                                acc=record.accuracy))
    else:
        fp.write(
            '>{name}\n{seq}\n'.format(name=record.name,
                                      seq=record.sequence))


def write_record_pair(read1, read2, fp):
    """
    Writes pair of sequence records to 'fp' in FASTA/FASTQ format.
    """
    if hasattr(read1, 'accuracy'):
        assert hasattr(read2, 'accuracy')
    write_record(read1, fp)
    write_record(read2, fp)


# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
