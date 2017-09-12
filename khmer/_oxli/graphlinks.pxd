from libcpp.memory cimport shared_ptr

from khmer._oxli.oxli_types cimport *
from khmer._oxli.graphs cimport CpHashgraph, Hashgraph, Nodegraph, Countgraph


cdef class Link:
    cdef readonly HashIntoType u
    cdef readonly HashIntoType v

    @staticmethod
    cdef Link _new(HashIntoType u, HashIntoType v)


cdef class LinkPath:

    cdef readonly list path


cdef class GraphLinker:

    cdef Hashgraph graph
    cdef shared_ptr[CpHashgraph] _graph
    cdef dict links

    cdef list _get_junction_choices(self, string sequence)
