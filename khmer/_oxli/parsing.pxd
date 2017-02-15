from libcpp.memory cimport unique_ptr
from libcpp cimport bool
from libcpp.string cimport string

from wrapper cimport *


cdef class Alphabets:

    @staticmethod
    cdef string _get(string name)


cdef class Sequence:
    cdef CpSequence _obj

    @staticmethod
    cdef Sequence _wrap(CpSequence cseq)


cdef class ReadBundle:
    cdef list reads


cdef class FastxParser:
    cdef unique_ptr[CpFastxParser] _this

    cpdef bool is_complete(self)
    cdef Sequence _next(self)


cdef class SanitizedFastxParser(FastxParser):
    cdef readonly int n_bad
    cdef readonly string _alphabet
    cdef bool convert_n

    cpdef bool is_complete(self)
    cdef Sequence _next(self)


cdef class SplitPairedReader:

    cdef FastxParser left_parser
    cdef FastxParser right_parser
    cdef readonly int min_length
    cdef readonly bool force_name_match

    cdef tuple _next(self)


cdef class BrokenPairedReader:

    cdef FastxParser parser
    cdef readonly int min_length
    cdef readonly bool force_single
    cdef readonly bool require_paired
    cdef readonly Sequence record

    cdef tuple _next(self)


cpdef tuple _split_left_right(str name)

cdef int _check_is_pair(Sequence first, Sequence second)

cpdef bool check_is_left(str name)

cpdef bool check_is_right(str name)

cdef inline bool is_valid(const char base, string& alphabet)

cdef inline bool sanitize_sequence(string& sequence,
                                   string& alphabet,
                                   bool convert_n)

