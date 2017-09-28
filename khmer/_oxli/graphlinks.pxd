cimport cython
from libcpp.memory cimport shared_ptr
from libcpp.list cimport list as stdlist
from libc.stdint cimport uint8_t, uint32_t, uint64_t

from khmer._oxli.oxli_types cimport *
from khmer._oxli.graphs cimport CpHashgraph, Hashgraph, Nodegraph, Countgraph


cdef extern from "oxli/links.hh" namespace "oxli":
    ctypedef struct Junction:
        HashIntoType u
        HashIntoType v
        uint64_t count

    ctypedef stdlist[Junction*] JunctionList

    cdef cppclass CpLink "oxli::Link":
        CpLink(uint64_t, bool)
        CpLink(uint64_t)
        bool is_forward()
        Junction* start_junction()
        Junction* end_junction()

        stdlist[Junction*].iterator begin()
        stdlist[Junction*].iterator end()
        const JunctionList& get_junctions()

    ctypedef stdlist[CpLink*] LinkList
    
    cdef cppclass CpGraphLinker "oxli::GraphLinker":
        CpGraphLinker(shared_ptr[CpHashgraph])

        Junction* get_junction(HashIntoType)
        Junction* get_junction(HashIntoType, HashIntoType)
        Junction* get_junction(Junction&)
        shared_ptr[JunctionList] get_junctions(const string&)

        void add_links(const string&)
        shared_ptr[LinkList] get_links(const string&)
        void report() const


cdef class GraphLinker:

    cdef shared_ptr[CpHashgraph] _graph
    cdef shared_ptr[CpGraphLinker] _gl_this
