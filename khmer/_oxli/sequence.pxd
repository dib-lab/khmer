from libcpp cimport bool
from libcpp.memory cimport shared_ptr
from libcpp.utility cimport pair
from libcpp.string cimport string



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
        void set_clean_seq()

    ctypedef pair[CpSequence,CpSequence] CpSequencePair \
        "oxli::read_parsers::ReadPair"


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
    cdef string _get(str name) except *


cdef class Sequence:
    cdef CpSequence _obj

    @staticmethod
    cdef Sequence _wrap(CpSequence cseq)


cdef class ReadBundle:
    cdef list reads

cdef bool is_valid(const char base, string& alphabet)

cdef bool sanitize_sequence(string& sequence,
                                   string& alphabet,
                                   bool convert_n)
