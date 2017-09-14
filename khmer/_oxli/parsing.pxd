# -*- coding: UTF-8 -*-


from libc.stdint cimport uintptr_t

from libcpp cimport bool
from libcpp.memory cimport unique_ptr, shared_ptr, weak_ptr
from libcpp.utility cimport pair
from libcpp.string cimport string

from khmer._oxli.utils cimport oxli_raise_py_error


'''
extern declarations for liboxli.
'''

# C++ ostream wrapper code stolen shamelessly from stackoverflow
# http://stackoverflow.com/questions/30984078/cython-working-with-c-streams
# We need ostream to wrap ReadParser
cdef extern from "<iostream>" namespace "std":
    cdef cppclass ostream:
        ostream& write(const char*, int) except +

# obviously std::ios_base isn't a namespace, but this lets
# Cython generate the connect C++ code
cdef extern from "<iostream>" namespace "std::ios_base":
    cdef cppclass open_mode:
        pass
    cdef open_mode binary
    # you can define other constants as needed


cdef extern from "<fstream>" namespace "std":
    cdef cppclass ofstream(ostream):
        # constructors
        ofstream(const char*) except +
        ofstream(const char*, open_mode) except+


cdef extern from  "oxli/read_parsers.hh" namespace "oxli::read_parsers":
    cdef cppclass CpSequence "oxli::read_parsers::Read":
        string name
        string description
        string sequence
        string quality
        string cleaned_seq

        void reset()
        void write_fastx(ostream&)
        void set_cleaned_seq()

    ctypedef pair[CpSequence,CpSequence] CpSequencePair \
        "oxli::read_parsers::ReadPair"

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


cdef extern from "oxli/alphabets.hh" namespace "oxli":
    cdef string DNA_SIMPLE "oxli::alphabets::DNA_SIMPLE"
    cdef string DNAN_SIMPLE "oxli::alphabets::DNAN_SIMPLE"
    cdef string RNA_SIMPLE "oxli::alphabets::RNA_SIMPLE"
    cdef string RNAN_SIMPLE "oxli::alphabets::RNAN_SIMPLE"
    cdef string IUPAC_NUCL "oxli::alphabets::IUPAC_NUCL"
    cdef string IUPAC_AA "oxli::alphabets::IUPAC_AA"

'''
Extension Classes wrapping liboxli.
'''

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

cdef inline bool is_valid(const char base, string& alphabet)

cdef inline bool sanitize_sequence(string& sequence,
                                   string& alphabet,
                                   bool convert_n)
