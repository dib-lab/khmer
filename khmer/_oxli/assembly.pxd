from libcpp.memory cimport unique_ptr
from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdint cimport uint16_t

from oxli_types cimport *
from hashing cimport CpKmer, Kmer
from graphs cimport (CpHashgraph, CpHashtable, CpLabelHash, 
                     get_hashgraph_ptr, get_labelhash_ptr)


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
    cdef unique_ptr[CpLinearAssembler] _this

    cdef public object graph
    cdef CpHashgraph * _graph_ptr

    cdef public object stop_filter
    cdef CpHashgraph * _stop_filter_ptr
    
    cdef str _assemble(self, Kmer start)
    cdef str _assemble_left(self, Kmer start)
    cdef str _assemble_right(self, Kmer start)


cdef class SimpleLabeledAssembler:
    cdef unique_ptr[CpSimpleLabeledAssembler] _this

    cdef public object labels
    cdef CpLabelHash * _label_ptr

    cdef public object stop_filter
    cdef CpHashgraph * _stop_filter_ptr
    
    cdef vector[string] _assemble(self, Kmer start)

'''
cdef class JunctionCountAssembler:
    cdef unique_ptr[CpJunctionCountAssembler] _this

    cdef public object graph
    cdef CpHashgraph * _graph_ptr

    cdef public object stop_filter
    cdef CpHashgraph * _stop_filter_ptr
    
    cdef str _assemble(self, Kmer)
    cdef uint16_t _consume(self, string)
'''
