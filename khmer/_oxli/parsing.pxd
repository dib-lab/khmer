# -*- coding: UTF-8 -*-


from libc.stdint cimport uintptr_t

from libcpp cimport bool
from libcpp.memory cimport unique_ptr, shared_ptr, weak_ptr
from libcpp.utility cimport pair
from libcpp.string cimport string

from khmer._oxli.utils cimport oxli_raise_py_error
from khmer._oxli.sequence cimport Sequence, CpSequence, CpSequencePair


'''
extern declarations for liboxli.
'''

cdef extern from  "oxli/read_parsers.hh" namespace "oxli::read_parsers":

    cdef cppclass CpReadParser "oxli::read_parsers::ReadParser" [SeqIO]:
        CpReadParser(unique_ptr[SeqIO]) except+
        CpReadParser(CpReadParser&)
        CpReadParser& operator=(CpReadParser&)
        CpReadParser(CpReadParser&&)
        CpReadParser& operator=(CpReadParser&&)

        CpSequence get_next_read()
        CpSequencePair get_next_read_pair()
        CpSequencePair get_next_read_pair(uint8_t)

        uintptr_t get_num_reads()
        bool is_complete()
        void close()

    cdef cppclass CpFastxReader "oxli::read_parsers::FastxReader":
        CpFastxReader() except+
        CpFastxReader(const string&) except+

        CpFastxReader(CpFastxReader&)
        CpFastxReader& operator=(CpFastxReader&)

        CpFastxReader(CpFastxReader&&)
        CpFastxReader& operator=(CpFastxReader&&)

        CpSequence get_next_read()
        bool is_complete()
        uintptr_t get_num_reads()
        void close()


    shared_ptr[CpReadParser[SeqIO]] get_parser[SeqIO](const string&) except +oxli_raise_py_error
    ctypedef shared_ptr[CpReadParser[CpFastxReader]] FastxParserPtr
    ctypedef weak_ptr[CpReadParser[CpFastxReader]] WeakFastxParserPtr


cdef extern from "khmer/_cpy_khmer.hh":
    ctypedef struct CPyReadParser_Object "khmer::khmer_ReadParser_Object":
        FastxParserPtr parser



cdef class FastxParser:
    cdef shared_ptr[CpReadParser[CpFastxReader]] _this

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


cpdef tuple _split_left_right(unicode s)

cdef tuple _cppstring_split_left_right(string& s)

cdef int _check_is_pair(Sequence first, Sequence second)

cpdef bool check_is_left(s)

cpdef bool check_is_right(s)
