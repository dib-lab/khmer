from libcpp.memory cimport unique_ptr
from libcpp cimport bool

from _oxli cimport CpSequence, CpFastxParser

cdef class Sequence:
    cdef unique_ptr[CpSequence] _this


cdef class SequencePair:
    cdef Sequence first
    cdef Sequence second


cdef class ReadBundle:
    cdef list reads


cdef class FastxParser:
    cdef unique_ptr[CpFastxParser] _this
    cdef Sequence _next(self)


cpdef tuple broken_paired_reader(FastxParser fastx_iter, int min_length=*,
                                 bool force_single=*, bool require_paired=*)

cdef tuple _split_left_right(str name)

cdef bool check_is_pair(Sequence first, Sequence second)

cdef bool check_is_left(str name)

cdef bool check_is_right(str name)


