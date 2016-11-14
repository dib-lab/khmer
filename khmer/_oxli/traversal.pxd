from ._oxli cimport *

cdef extern from "traversal.hh":
    cdef cppclass Traverser "khmer::Traverser":
        Traverser(Hashtable *)

        void push_filter(KmerFilter)
        KmerFilter pop_filter()
    
        uint32_t traverse(const Kmer&, KmerQueue&) const
        uint32_t traverse_left(const Kmer&, KmerQueue&) const     
        uint32_t traverse_right(const Kmer&, KmerQueue&) const

        uint32_t degree(const Kmer&) const
        uint32_t degree_left(const Kmer&) const
        uint32_t degree_right(const Kmer&) const



