from cython.operator cimport dereference as deref
from libcpp.memory cimport make_shared

from khmer._oxli.utils cimport _bstring

cdef class GraphLinker:

    def __cinit__(self, Hashgraph graph):
        self._graph = graph._hg_this
        
        if type(self) is GraphLinker:
            self._gl_this = make_shared[CpGraphLinker](self._graph)

    def add_links(self, str sequence):
        cdef string _sequence = _bstring(sequence)
        deref(self._gl_this).add_links(_sequence)

    def get_junctions(self, str sequence):
        cdef string _sequence = _bstring(sequence)
        cdef shared_ptr[JunctionList] junctions = deref(self._gl_this).get_junctions(_sequence)
        cdef Junction* j
        for j in deref(junctions):
            yield deref(j)

    def get_links(self, str sequence):
        cdef string _sequence = _bstring(sequence)
        cdef shared_ptr[LinkList] links = deref(self._gl_this).get_links(_sequence)
        
        cdef CpLink* link
        cdef stdlist[Junction*].iterator it
        for link in deref(links):
            ret = []
            it = deref(link).begin()
            while it != deref(link).end():
                ret.append(deref(deref(it)))
            yield ret
        

    def report(self):
        deref(self._gl_this).report()
