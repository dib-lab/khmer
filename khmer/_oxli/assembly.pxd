from libcpp.memory cimport shared_ptr, make_shared
from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdint cimport uint16_t

from khmer._oxli.oxli_types cimport *
from khmer._oxli.hashing cimport CpKmer, Kmer
from khmer._oxli.graphs cimport Hashgraph, CpHashgraph, CpHashtable
from khmer._oxli.labeling cimport CpLabelHash, GraphLabels


cdef extern from "oxli/assembler.hh" namespace "oxli":
    cdef cppclass CpLinearAssembler "oxli::LinearAssembler":
        CpLinearAssembler(CpHashgraph *)
    
        string assemble(const CpKmer, const CpHashgraph *) const
        string assemble_left(const CpKmer, const CpHashgraph *) const     
        string assemble_right(const CpKmer, const CpHashgraph *) const

        string assemble(const CpKmer) const
        string assemble_left(const CpKmer) const     
        string assemble_right(const CpKmer) const

    cdef cppclass CpSimpleLabeledAssembler "oxli::SimpleLabeledAssembler":
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


cdef class LinearAssembler:
    cdef shared_ptr[CpLinearAssembler] _this

    cdef public Hashgraph graph
    cdef shared_ptr[CpHashgraph] _graph_ptr

    cdef public Hashgraph stop_filter
    cdef shared_ptr[CpHashgraph] _stop_filter_ptr
    
    cdef str _assemble(self, CpKmer start)
    cdef str _assemble_left(self, CpKmer start)
    cdef str _assemble_right(self, CpKmer start)


cdef class SimpleLabeledAssembler:
    cdef shared_ptr[CpSimpleLabeledAssembler] _this

    cdef public GraphLabels labels
    cdef public Hashgraph graph
    cdef shared_ptr[CpLabelHash] _label_ptr

    cdef public Hashgraph stop_filter
    cdef shared_ptr[CpHashgraph] _stop_filter_ptr
    
    cdef vector[string] _assemble(self, CpKmer start)


cdef class JunctionCountAssembler:
    cdef shared_ptr[CpJunctionCountAssembler] _this

    cdef public Hashgraph graph
    cdef shared_ptr[CpHashgraph] _graph_ptr

    cdef public Hashgraph stop_filter
    cdef shared_ptr[CpHashgraph] _stop_filter_ptr
    
    cdef vector[string] _assemble(self, CpKmer)
