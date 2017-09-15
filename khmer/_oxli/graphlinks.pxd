cimport cython
from libcpp.memory cimport shared_ptr
from libc.stdint cimport uint8_t, uint32_t, uint64_t

from khmer._oxli.oxli_types cimport *
from khmer._oxli.graphs cimport CpHashgraph, Hashgraph, Nodegraph, Countgraph


cdef enum DBGNucl:
    A, C, G, T


cdef class Junction:
    cdef readonly HashIntoType u
    cdef readonly HashIntoType v

    @staticmethod
    cdef Junction _new(HashIntoType u, 
                   HashIntoType v)


cdef class Link:

    cdef readonly list junctions
    cdef bool forward

    @staticmethod
    cdef Link _new(list junctions, bool forward)


cdef class GraphLinker:

    cdef Hashgraph graph
    cdef shared_ptr[CpHashgraph] _graph
    cdef WordLength K
    cdef dict links

    cdef tuple _get_junctions(self, string sequence)
    cdef int _add_link(self, string sequence,
            unsigned int min_link_size=*) except -1
