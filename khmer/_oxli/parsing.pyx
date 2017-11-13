# -*- coding: UTF-8 -*-


from cython.operator cimport dereference as deref
cimport cython
from libcpp cimport bool
from libcpp.string cimport string

import sys

from khmer._oxli.utils cimport _bstring, _ustring


cdef class Alphabets:
    
    @staticmethod
    def get(name):
        cdef unicode alphabet = _ustring(Alphabets._get(_bstring(name)))
        if not alphabet:
            raise ValueError('No alphabet with name {0}'.format(name))
        return alphabet

    @staticmethod
    cdef string _get(string name):
        if name == b'DNA_SIMPLE':
            return DNA_SIMPLE
        elif name == b'DNAN_SIMPLE':
            return DNAN_SIMPLE
        elif name == b'RNA_SIMPLE':
            return RNA_SIMPLE
        elif name == b'RNAN_SIMPLE':
            return RNAN_SIMPLE
        elif name == b'IUPAC_NUCL':
            return IUPAC_NUCL
        elif name == b'IUPAC_AA':
            return IUPAC_AA
        else:
            return string()


@cython.freelist(100)
cdef class Sequence:

    def __cinit__(self, name=None, sequence=None,
                        quality=None, description=None,
                        cleaned_seq=None):

        if name is not None and sequence is not None:
            self._obj.sequence = _bstring(sequence)
            self._obj.name = _bstring(name)
            if description is not None:
                self._obj.description = _bstring(description)
            if quality is not None:
                self._obj.quality = _bstring(quality)
            if cleaned_seq is not None:
                self._obj.cleaned_seq = _bstring(cleaned_seq)
            else:
                self._obj.cleaned_seq = self._obj.sequence

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return 'Sequence(name="{0}", sequence="{1}")'.format(self.name, self.sequence)

    def __len__(self):
        return self._obj.sequence.length()

    def __richcmp__(x, y, op):
        if op == 2:
            return x.name == y.name and x.sequence == y.sequence
        else:
            raise NotImplementedError('Operator not available')

    def kmers(self, int K):
        cdef int i = 0
        cdef unicode sequence = self.sequence
        for i in range(0, len(self)-K+1):
            yield sequence[i:i+K]

    def __getitem__(self, x):
        # Definitely optimize this.
        return self.sequence[x]

    @property
    def name(self):
        cdef unicode name = self._obj.name
        return self._obj.name if name else None

    @property
    def sequence(self):
        cdef unicode sequence = self._obj.sequence
        return self._obj.sequence if sequence else None

    @property
    def description(self):
        cdef unicode description = self._obj.description
        return description if description else None

    @property
    def quality(self):
        cdef unicode quality = self._obj.quality
        return quality if quality else None

    @property
    def cleaned_seq(self):
        cdef unicode cleaned_seq = self._obj.cleaned_seq
        return cleaned_seq if cleaned_seq else None

    @staticmethod
    def from_screed_record(record):
        cdef Sequence seq = Sequence(name=record.name,
                                     sequence=record.sequence)
        if hasattr(record, 'quality'):
            seq._obj.quality = _bstring(record.quality)

        for attr in ('annotations', 'description'):
            if hasattr(record, attr):
                seq._obj.description = _bstring(getattr(record, attr))

        return seq

    @staticmethod
    cdef Sequence _wrap(CpSequence cseq):
        cdef Sequence seq = Sequence()
        seq._obj = cseq
        return seq


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


cdef inline bool is_valid(const char base, string& alphabet):
    cdef char b
    for b in alphabet:
        if b == base:
            return True
    return False


cdef inline bool sanitize_sequence(string& sequence,
                                   string& alphabet,
                                   bool convert_n):
    cdef int i = 0
    for i in range(sequence.length()):
        sequence[i] &= 0xdf
        if not is_valid(sequence[i], alphabet):
            return False
        if convert_n and sequence[i] == b'N':
            sequence[i] = b'A'
    return True


cdef class FastxParser:

    def __cinit__(self, filename, *args, **kwargs):
        self._this = get_parser[CpFastxReader](_bstring(filename))

    cdef Sequence _next(self):
        if not self.is_complete():
            return Sequence._wrap(deref(self._this).get_next_read())
        else:
            return None

    cpdef bool is_complete(self):
        return deref(self._this).is_complete()

    def __iter__(self):
        cdef Sequence seq
        while not self.is_complete():
            seq = self._next()
            yield seq


cdef class SanitizedFastxParser(FastxParser):

    def __cinit__(self, filename, alphabet='DNAN_SIMPLE',
                        bool convert_n=True):
        self.n_bad = 0
        self.convert_n = convert_n
        self._alphabet = Alphabets._get(_bstring(alphabet))

    cdef Sequence _next(self):
        cdef Sequence seq
        cdef bool good

        if not self.is_complete():
            seq = FastxParser._next(self)
            good = sanitize_sequence(seq._obj.sequence,
                                     self._alphabet,
                                     self.convert_n)
            if not good:
                self.n_bad += 1
                return None
            else:
                return seq
        else:
            return None

    def __iter__(self):
        cdef Sequence seq
        while not self.is_complete():
            seq = self._next()
            if seq is not None:
                yield seq


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
                if len(first) >= self.min_length and \
                   len(second) >= self.min_length:

                    yield read_num, True, first, second
            else:
                yield read_num, True, first, second

            read_num += found
            found, first, second, err = self._next()

    cdef tuple _next(self):
        cdef Sequence first = self.left_parser._next()
        cdef bool first_complete = self.left_parser.is_complete()

        cdef Sequence second = self.right_parser._next()
        cdef bool second_complete = self.right_parser.is_complete()
        

        if first_complete is not second_complete:
            err = UnpairedReadsError('Differing lengths of left '\
                                     'and right files!')
            return -1, None, None, err

        if first_complete:
            return 0, None, None, None

        if first is None or second is None:
            return 1, first, second, None

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

        found, first, second, err = self._next()
        while (found != 0):
            if err is not None:
                raise err

            if self.min_length > 0:
                if first is not None and len(first) < self.min_length:
                    first = None
                    found -= 1
                if second is not None and len(second) < self.min_length:
                    second = None
                    found -= 1

            if self.force_single:
                if first is not None:
                    yield read_num, found == 2, first, None
                    read_num += found
                if second is not None:
                    yield read_num, found == 2, second, None
                    read_num += found
            elif self.require_paired:
                if first is not None and second is not None:
                    yield read_num, found == 2, first, second
                    read_num += found
            else:
                if first is not None or second is not None:
                    yield read_num, found == 2, first, second
                    read_num += found
            found, first, second, err = self._next()

    cdef tuple _next(self):
        cdef Sequence first, second
        cdef int is_pair

        if self.record is None:
            first = self.parser._next()
            if first is None:
                if self.parser.is_complete():
                    return 0, None, None, None
                else:
                    if self.require_paired:
                        err = UnpairedReadsError(
                            "Uneven number of reads when require_paired is set!",
                            first)
                        return -1, None, None, err
                    return 1, first, None, None
        else:
            first = self.record
        
        second = self.parser._next()
        
        # check if paired
        if second is not None and first is not None:
            is_pair = _check_is_pair(first, second)
            if is_pair == -1:
                err = ValueError("records must be same type (FASTA or FASTQ)")
                return -1, None, None, err
            if is_pair and not self.force_single:
                self.record = None
                return 2, first, second, None    # found 2 proper records
            else:   # orphan.
                if self.require_paired:
                    err = UnpairedReadsError(
                        "Unpaired reads when require_paired is set!",
                        first, second)
                    return -1, None, None, err
                self.record = second
                return 1, first, None, None
        elif self.parser.is_complete():
            # ran out of reads getting second, handle last record
            if self.require_paired:
                err =  UnpairedReadsError("Unpaired reads when require_paired "
                                          "is set!", first, None)
                return -1, None, None, err
            self.record = None
            return 1, first, second, None
        else: # one read was invalid, but that doesn't mean they were unpaired
            return 1, first, second, None


cpdef tuple _split_left_right(unicode s):
    cdef string cppstr = s.encode('UTF-8')
    return _cppstring_split_left_right(cppstr)


cdef tuple _cppstring_split_left_right(string& s):
    """Split record name at the first whitespace and return both parts.

    RHS is set to an empty string if not present.
    """
    cdef unsigned int i
    cdef const char * c_str = s.c_str()
    cdef unicode lhs, rhs
    lhs = u''
    rhs = u''
    for i in range(len(s)):
        if lhs == u'':
            if (c_str[i] == b' ' or c_str[i] == b'\t'):
                lhs = _ustring(c_str[0:i])
        else:
            if c_str[i] != b' ' or c_str[i] != b'\t':
                rhs = _ustring(c_str[i:len(s)])
                break
    lhs = _ustring(c_str[0:len(s)])  if lhs == u'' else lhs
    return lhs, rhs


cdef int _check_is_pair(Sequence first, Sequence second):
    """Check if the two sequence records belong to the same fragment.

    In an matching pair the records are left and right pairs
    of each other, respectively.  Returns True or False as appropriate.

    Handles both Casava formats: seq/1 and seq/2, and 'seq::... 1::...'
    and 'seq::... 2::...'.

    Also handles the default format of the SRA toolkit's fastq-dump:
    'Accession seq/1'
    """
    if first.quality is None or second.quality is None:
        if first.quality is not second.quality:
            return -1

    cdef unicode lhs1, rhs1, lhs2, rhs2
    lhs1, rhs1 = _cppstring_split_left_right(first._obj.name)
    lhs2, rhs2 = _cppstring_split_left_right(second._obj.name)

    # handle 'name/1'
    cdef unicode subpart1, subpart2
    if lhs1.endswith('/1') and lhs2.endswith('/2'):
        subpart1 = lhs1.split('/', 1)[0]
        subpart2 = lhs2.split('/', 1)[0]

        if subpart1 and subpart1 == subpart2:
            return 1

    # handle '@name 1:rst'
    elif lhs1 == lhs2 and rhs1.startswith('1:') and rhs2.startswith('2:'):
        return 1

    # handle @name seq/1
    elif lhs1 == lhs2 and rhs1.endswith('/1') and rhs2.endswith('/2'):
        subpart1 = rhs1.split('/', 1)[0]
        subpart2 = rhs2.split('/', 1)[0]

        if subpart1 and subpart1 == subpart2:
            return 1

    return 0


def check_is_pair(first, second):
    if type(first) is not Sequence:
        first = Sequence.from_screed_record(first)
    if type(second) is not Sequence:
        second = Sequence.from_screed_record(second)
    cdef int ret = _check_is_pair(first, second)
    if ret == -1:
        raise ValueError("both records must be same type (FASTA or FASTQ)")
    return ret == 1


cpdef bool check_is_left(s):
    """Check if the name belongs to a 'left' sequence (/1).

    Returns True or False.

    Handles both Casava formats: seq/1 and 'seq::... 1::...'
    """
    cdef unicode lhs, rhs
    lhs, rhs = _split_left_right(_ustring(s))
    if lhs.endswith('/1'):              # handle 'name/1'
        return True
    elif rhs.startswith('1:'):          # handle '@name 1:rst'
        return True

    elif rhs.endswith('/1'):            # handles '@name seq/1'
        return True

    return False


cpdef bool check_is_right(s):
    """Check if the name belongs to a 'right' sequence (/2).

    Returns True or False.

    Handles both Casava formats: seq/2 and 'seq::... 2::...'
    """
    cdef unicode lhs, rhs
    lhs, rhs = _split_left_right(_ustring(s))
    if lhs.endswith('/2'):              # handle 'name/2'
        return True
    elif rhs.startswith('2:'):          # handle '@name 2:rst'
        return True

    elif rhs.endswith('/2'):            # handles '@name seq/2'
        return True

    return False

