from libcpp.memory cimport shared_ptr
from wrapper cimport CpKmer, HashIntoType, WordLength

cdef class Kmer:
    cdef shared_ptr[CpKmer] _this
    cdef readonly str kmer

    @staticmethod
    cdef Kmer wrap(CpKmer * cpkmer, WordLength K)
