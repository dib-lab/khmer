# cython: c_string_type=unicode, c_string_encoding=utf8
from __future__ import print_function
from cython.operator cimport dereference as deref
from libcpp cimport bool
from libcpp.string cimport string

import sys

from _oxli cimport CpSequence


cdef inline void move_cpsequence(CpSequence& move_from, CpSequence& move_to):
    move_to.name.assign(move_from.name)
    move_to.sequence.assign(move_from.sequence)
    move_to.annotations.assign(move_from.annotations)
    move_to.quality.assign(move_from.quality)

cdef inline void move_sequence(Sequence move_from, Sequence move_to):
    move_cpsequence(deref(move_from._this), deref(move_to._this))

cdef class Sequence:

    def __cinit__(self, create=False):
        if create is True:
            self._this.reset(new CpSequence())

    cdef void take(self, Sequence other):
        move_cpsequence(deref(other._this), deref(self._this))

    def __str__(self):
        return self.sequence

    def __len__(self):
        return deref(self._this).sequence.length()

    @property
    def name(self):
        return deref(self._this).name

    @property
    def sequence(self):
        return deref(self._this).sequence

    @property
    def annotations(self):
        cdef str annotations = deref(self._this).annotations
        return annotations if annotations else None

    @property
    def quality(self):
        cdef str quality = deref(self._this).quality
        return quality if quality else None

    @staticmethod
    def create(str name, str sequence, str annotations=None, str quality=None):
        cdef Sequence seq = Sequence(create=True)
        deref(seq._this).sequence = sequence.encode('UTF-8')
        deref(seq._this).name = name.encode('UTF-8')
        if annotations is not None:
            deref(seq._this).annotations = annotations.encode('UTF-8')
        if quality is not None:
            deref(seq._this).quality = quality.encode('UTF-8')

        return seq

    @staticmethod
    cdef Sequence _create(str name, str sequence, str annotations=None, str quality=None):
        cdef Sequence seq = Sequence(create=True)
        deref(seq._this).sequence = sequence.encode('UTF-8')
        deref(seq._this).name = name.encode('UTF-8')
        if annotations is not None:
            deref(seq._this).annotations = annotations.encode('UTF-8')
        if quality is not None:
            deref(seq._this).quality = quality.encode('UTF-8')

        return seq


cdef class SequencePair:
    
    def __cinit__(self, Sequence first, Sequence second):
        self.first = first
        self.second = second


cdef class ReadBundle:

    def __cinit__(self, *raw_records):
        self.reads = [r for r in raw_records if r]

    @property
    def num_reads(self):
        return len(self.reads)

    @property
    def total_length(self):
        return sum([len(r.sequence) for r in self.reads])


def print_error(msg):
    """Print the given message to 'stderr'."""

    print(msg, file=sys.stderr)


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


cdef inline bool is_valid_dna(const char base):
    return base == 'A' or base == 'C' or base == 'G'\
           or base == 'T' or base == 'N'


cdef inline bool sanitize_sequence(string& sequence):
    cdef int i = 0
    for i in range(sequence.length()):
        sequence[i] &= 0xdf
        if not is_valid_dna(sequence[i]):
            return False
        if sequence[i] == 'N':
            sequence[i] = 'A'
    return True


cdef class FastxParser:

    def __cinit__(self, str filename):
        self._this.reset(new CpFastxParser(filename.encode()))

    cdef bool _cp_imprint_next(self, CpSequence& seq):
        if deref(self._this).is_complete():
            return False
        else:
            deref(self._this).imprint_next_read(seq)
            return True

    cdef bool _imprint_next(self, Sequence seq):
        if deref(self._this).is_complete():
            return False
        else:
            deref(self._this).imprint_next_read(deref(seq._this))
            return True

    cdef Sequence _next(self):
        cdef Sequence seq = Sequence(create=True)
        cdef bool res = self._imprint(seq)
        if res:
            return seq
        else:
            return None

    def __iter__(self):
        cdef Sequence seq = self._next()
        while seq is not None:
            yield seq
            seq = self._next()


cdef class SplitPairedReader:

    def __cinit__(self, FastxParser left_parser,
                         FastxParser right_parser,
                         int min_length=-1,
                         bool force_name_match=False):

        self.left_parser = left_parser
        self.right_parser = right_parser
        self.min_length = min_length
        self.force_name_match = force_name_match

    def __iter__(self):
        cdef Sequence first, second
        cdef object err
        cdef read_num = 0
        cdef int found

        found, first, second, err = self._next()
        while found != 0:
            if err is not None:
                raise err
            
            if self.min_length > 0:
                if len(first) >= self.min_length or \
                   len(second) >= self.min_length:

                    yield read_num, True, first, second
            else:
                yield read_num, True, first, second

            read_num += 2
            found, first, second, err = self._next()

    cdef tuple _next(self):
        cdef Sequence first = Sequence(create=True)
        cdef Sequence second = Sequence(create=True)
        cdef bool left_has_next, right_has_next

        left_has_next = self.left_parser._imprint_next(first)
        right_has_next = self.right_parser._imprint_next(second)

        if left_has_next != right_has_next:
            err = UnpairedReadsError('Differing lengths of left '\
                                     'and right files!')
            return -1, None, None, err

        if left_has_next == False:
            return 0, None, None, None

        if self.force_name_match:
            if _check_is_pair(first, second):
                return 2, first, second, None
            else:
                err =  UnpairedReadsError('', first, second)
                return -1, None, None, err
        else:
            return 2, first, second, None

cdef class BrokenPairedReader:

    def __cinit__(self, FastxParser parser, 
                  int min_length=-1,
                  bool force_single=False, 
                  bool require_paired=False):
        
        if force_single and require_paired:
            raise ValueError("force_single and require_paired cannot both be set!")

        self.parser = parser
        self.min_length = min_length
        self.force_single = force_single
        self.require_paired = require_paired

        self.record = None

    def __iter__(self):
        cdef Sequence first
        cdef Sequence second
        cdef object err
        cdef int found
        cdef int read_num = 0
        cdef bool passed_length = True

        found, first, second, err = self._next()
        while (found != 0):
            if err is not None:
                raise err

            if self.min_length > 0:
                if found == 1:
                    if len(first) < self.min_length:
                        passed_length = False
                else:
                    if len(first) < self.min_length or len(second) < self.min_length:
                        passed_length = False

                if passed_length:
                    yield read_num, found == 2, first, second
                    read_num += found
            else:
                yield read_num, found == 2, first, second
                read_num += found
            found, first, second, err = self._next()
            passed_length = True

    cdef tuple _next(self):
        cdef bool has_next
        cdef Sequence first, second

        if self.record == None:
            # No previous sequence. Try to get one.
            first = Sequence(create=True)
            has_next = self.parser._imprint_next(first)
            # And none left? We're outta here.
            if not has_next:
                return 0, None, None, None
        else:
            first = self.record
        
        second = Sequence(create=True)
        has_next = self.parser._imprint_next(second)
        
        # check if paired
        if has_next:
            if _check_is_pair(first, second) and not self.force_single:
                self.record = None
                return 2, first, second, None    # found 2 proper records
            else:                                   # orphan.
                if self.require_paired:
                    err = UnpairedReadsError(
                        "Unpaired reads when require_paired is set!",
                        first, second)
                    return -1, None, None, err
                self.record = second
                return 1, first, None, None
        else:
            # ran out of reads, handle last record
            if self.require_paired:
                err =  UnpairedReadsError("Unpaired reads when require_paired "
                                          "is set!", first, None)
                return -1, None, None, err
            self.record = None
            return 1, first, None, None


cdef tuple _split_left_right(str name):
    """Split record name at the first whitespace and return both parts.

    RHS is set to an empty string if not present.
    """
    cdef list parts = name.split(None, 1)
    cdef str lhs = parts[0]
    cdef str rhs = parts[1] if len(parts) > 1 else ''
    return lhs, rhs


cdef bool _check_is_pair(Sequence first, Sequence second):
    """Check if the two sequence records belong to the same fragment.

    In an matching pair the records are left and right pairs
    of each other, respectively.  Returns True or False as appropriate.

    Handles both Casava formats: seq/1 and seq/2, and 'seq::... 1::...'
    and 'seq::... 2::...'.

    Also handles the default format of the SRA toolkit's fastq-dump:
    'Accession seq/1'
    """
    #if hasattr(first, 'quality') or hasattr(second, 'quality'):
    #    if not (hasattr(first, 'quality') and hasattr(first, 'quality')):
    #        raise ValueError("both records must be same type (FASTA or FASTQ)")

    cdef str lhs1, rhs1, lhs2, rhs2
    lhs1, rhs1 = _split_left_right(first.name)
    lhs2, rhs2 = _split_left_right(second.name)

    # handle 'name/1'
    cdef str subpart1, subpart2
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


def check_is_pair(Sequence first, Sequence second):
    cdef int res = _check_is_pair(first, second)
    return res == 1


cdef bool check_is_left(str name):
    """Check if the name belongs to a 'left' sequence (/1).

    Returns True or False.

    Handles both Casava formats: seq/1 and 'seq::... 1::...'
    """
    cdef str lhs, rhs
    lhs, rhs = _split_left_right(name)
    if lhs.endswith('/1'):              # handle 'name/1'
        return True
    elif rhs.startswith('1:'):          # handle '@name 1:rst'
        return True

    elif rhs.endswith('/1'):            # handles '@name seq/1'
        return True

    return False


cdef bool check_is_right(str name):
    """Check if the name belongs to a 'right' sequence (/2).

    Returns True or False.

    Handles both Casava formats: seq/2 and 'seq::... 2::...'
    """
    cdef str lhs, rhs
    lhs, rhs = _split_left_right(name)
    if lhs.endswith('/2'):              # handle 'name/2'
        return True
    elif rhs.startswith('2:'):          # handle '@name 2:rst'
        return True

    elif rhs.endswith('/2'):            # handles '@name seq/2'
        return True

    return False

