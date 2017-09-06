from libcpp cimport bool
from libcpp.memory cimport unique_ptr
from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdint cimport uint64_t, uint8_t

from khmer._oxli.oxli_types cimport *
from khmer._oxli.parsing cimport CpReadParser
from khmer._oxli.utils cimport oxli_raise_py_error


cdef extern from "oxli/hllcounter.hh" namespace "oxli":
    cdef cppclass CpHLLCounter "oxli::HLLCounter":
        CpHLLCounter(double, WordLength) except +oxli_raise_py_error
        CpHLLCounter(int, WordLength) except +oxli_raise_py_error

        void add(const string &)
        unsigned int consume_string(const string &)
        void consume_seqfile[SeqIO](const string &,
                                    bool,
                                    unsigned int &,
                                    uint64_t &) except +oxli_raise_py_error

#        void consume_seqfile[SeqIO](unique_ptr[CpReadParser[SeqIO]]&,
#                                    bool,
#                                    unsigned int &,
#                                    unsigned long long &)
        unsigned int check_and_process_read(string &, bool &)
        bool check_and_normalize_read(string &) const
        uint64_t estimate_cardinality()
        void merge(CpHLLCounter &) except +oxli_raise_py_error
        double get_alpha()
        int get_p()
        int get_ncounters()
        void set_ksize(WordLength) except +oxli_raise_py_error
        int get_ksize()
        vector[uint8_t] get_counters()
        void set_counters(vector[uint8_t]) except +oxli_raise_py_error
        double get_erate()
        void set_erate(double) except +oxli_raise_py_error


cdef class HLLCounter:
    cdef unique_ptr[CpHLLCounter] _this
    cpdef tuple consume_seqfile(self, filename, bool stream_records=*)
