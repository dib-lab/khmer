from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.queue cimport queue
from libcpp.memory cimport unique_ptr, weak_ptr, shared_ptr
from libcpp.utility cimport pair
from libc.stdint cimport uint32_t, uint8_t, uint16_t, uint64_t, uintptr_t 

from oxli_types cimport *
from hashing cimport CpKmer, KmerQueue, KmerSet, KmerFilter
from graphs cimport CpHashtable, CpHashgraph, CpHashtable, CpLabelHash
from parsing cimport CpReadParser, CpSequence, CpFastxReader


cdef extern from "oxli/assembler.hh" namespace "oxli":
    cdef cppclass CpLinearAssembler "oxli::LinearAssembler":
        CpLinearAssembler(CpHashgraph *)
    
        string assemble(const CpKmer, const CpHashgraph *) const
        string assemble_left(const CpKmer, const CpHashgraph *) const     
        string assemble_right(const CpKmer, const CpHashgraph *) const

        string assemble(const CpKmer) const
        string assemble_left(const CpKmer) const     
        string assemble_right(const CpKmer) const

    cdef cppclass CpSimpleLabeledAsssembler "oxli::SimpleLabeledAsssembler":
        CpSimpleLabeledAssembler(const CpLabelHash *)

        vector[string] assemble(const CpKmer)
        vector[string] assemble(const CpKmer, const CpHashgraph *) const

    cdef cppclass CpJunctionCountAssembler "oxli::JunctionCountAssembler":
        CpJunctionCountAssembler(CpHashgraph *)

        vector[string] assemble(const CpKmer) const
        vector[string] assemble(const CpKmer, const CpHashtable *) const
        uint16_t consume(string)
        void count_junction(CpKmer, CpKmer)
        BoundedCounterType get_junction_count(CpKmer, CpKmer) const


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

