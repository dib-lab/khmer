from libc.stdint cimport uint32_t
from libcpp.memory cimport shared_ptr

from hashing cimport CpKmer, KmerFilter, KmerQueue
from graphs cimport CpHashgraph

cdef extern from "oxli/traversal.hh" namespace "oxli":
    cdef cppclass CpTraverser "oxli::Traverser":
        CpTraverser(CpHashgraph *)

        void push_filter(KmerFilter)
        KmerFilter pop_filter()
    
        uint32_t traverse(const CpKmer&, KmerQueue&) const
        uint32_t traverse_left(const CpKmer&, KmerQueue&) const     
        uint32_t traverse_right(const CpKmer&, KmerQueue&) const

        uint32_t degree(const CpKmer&) const
        uint32_t degree_left(const CpKmer&) const
        uint32_t degree_right(const CpKmer&) const


cdef class Traverser:
    cdef shared_ptr[CpTraverser] _this
    cdef shared_ptr[CpHashgraph] _graph_ptr

