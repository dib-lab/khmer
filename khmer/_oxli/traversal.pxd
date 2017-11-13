from libc.stdint cimport uint32_t
from libcpp.memory cimport shared_ptr
from libcpp cimport bool

from khmer._oxli.hashing cimport Kmer, CpKmer, KmerFilter, KmerQueue
from khmer._oxli.graphs cimport Hashgraph, CpHashgraph


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
    cdef Hashgraph graph
    cdef shared_ptr[CpTraverser] _this
    cdef shared_ptr[CpHashgraph] _graph_ptr
    cdef list _kmerqueue_to_kmer_list(self, KmerQueue * kmers)
    cdef list _kmerqueue_to_hash_list(self, KmerQueue * kmers)
    cdef list _neighbors(self, CpKmer start, int direction=*)
