from libcpp.memory cimport unique_ptr
from libcpp cimport bool
from libcpp.string cimport string

from _oxli cimport CpSequence, CpFastxParser

cdef inline void move_cpsequence(CpSequence& move_from, CpSequence& move_to)
cdef inline void move_sequence(Sequence move_from, Sequence move_to)

cdef class Sequence:
    cdef unique_ptr[CpSequence] _this

    cdef void take(self, Sequence other)
    @staticmethod
    cdef Sequence _create(str name, str sequence, 
                          str annotations=*, str quality=*)

cdef class SequencePair:
    cdef Sequence first
    cdef Sequence second


cdef class ReadBundle:
    cdef list reads


cdef class FastxParser:
    cdef unique_ptr[CpFastxParser] _this

    cdef bool _imprint_next(self, Sequence seq)
    cdef bool _cp_imprint_next(self, CpSequence& seq)
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

cdef tuple _split_left_right(str name)

cdef bool _check_is_pair(Sequence first, Sequence second)

cdef bool check_is_left(str name)

cdef bool check_is_right(str name)

cdef inline bool is_valid_dna(const char base)

cdef inline bool sanitize_sequence(string& sequence)

