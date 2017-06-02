from libcpp cimport bool
from libcpp.memory cimport unique_ptr
from libcpp.sting cimport string
from libcpp.vector cimport vector
from libc.stdint cimport uint64_t

from oxli_types cimport *
from parsing cimport CpReadParser

cdef extern from "oxli/hllcounter.hh" namespace "oxli":
    cdef cppclass CpHLLCounter "oxli::HLLCounter":
        CpHLLCounter(double, WordLength)
        CpHLLCounter(int, WordLength)

        void add(const string &)
        unsigned int consume_string(const string &)
        void consume_seqfile[SeqIO](const string &,
                                    bool,
                                    unsigned int &,
                                    unsigned long long &)

        void consume_seqfile[SeqIO](unique_ptr[CpReadParser[SeqIO]]&,
                                    bool,
                                    unsigned int &,
                                    unsigned long long &)
        unsigned int check_and_process_read(string &, bool &)
        bool check_and_normalize_read(string &) const
        uint64_t estimate_cardinality()
        void merge(CpHLLCounter &)
        double get_alpha()
        int get_p()
        int get_m()
        void set_ksize(WordLegth)
        int get_ksize()
        vector[int] get_M()
        double get_erate()
        void set_erate(double)

