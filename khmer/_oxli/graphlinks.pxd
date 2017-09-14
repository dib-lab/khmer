cimport cython
from libcpp.memory cimport shared_ptr
from libc.stdint cimport uint8_t, uint32_t, uint64_t

from khmer._oxli.oxli_types cimport *
from khmer._oxli.graphs cimport CpHashgraph, Hashgraph, Nodegraph, Countgraph


cdef enum DBGNucl:
    A, C, G, T


cdef class Link:
    cdef readonly HashIntoType u
    cdef readonly HashIntoType v
    cdef readonly bool forward
    cdef list children
    cdef Link parent

    cpdef bool add_child(self, DBGNucl nuc, Link child_link)
    cpdef Link get_child(self, DBGNucl nuc)

    @staticmethod
    cdef Link _new(HashIntoType u, 
                   HashIntoType v, 
                   bool forward,
                   Link parent=*)


cdef class LinkPath:

    cdef readonly list path


cdef class GraphLinker:

    cdef Hashgraph graph
    cdef shared_ptr[CpHashgraph] _graph
    cdef WordLength K
    cdef dict links

    cdef list _get_junction_choices(self, string sequence)
