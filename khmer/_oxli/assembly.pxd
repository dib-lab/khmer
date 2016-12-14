from libcpp.memory cimport unique_ptr
from libcpp.string import string

from _oxli cimport CpLinearAssembler, CpCompactingAssembler, CpHashtable, CpKmer
from hashing cimport Kmer

cdef class LinearAssembler:
    cdef unique_ptr[CpLinearAssembler] _this

    cdef public object graph
    cdef CpHashtable * _graph_ptr

    cdef public object stop_filter
    cdef CpHashtable * _stop_filter_ptr
    
    cdef str _assemble(self, Kmer start)
    cdef str _assemble_left(self, Kmer start)
    cdef str _assemble_right(self, Kmer start)


cdef class CompactingAssembler(LinearAssembler):
    pass
    #cdef unique_ptr[CpCompactingAssembler] _this

