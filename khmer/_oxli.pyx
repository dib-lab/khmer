import cython
from cython.operator cimport dereference as deref, preincrement as inc
from libc.limits cimport UINT_MAX
from libcpp.memory cimport unique_ptr, weak_ptr

from khmer._khmer import Countgraph
from khmer._khmer import Nodegraph

cdef class Component:

    cdef ComponentPtr _this

    def __cinit__(self, Component other=None):
        if other is not None:
            self._this.reset(other._this.get())
        else:
            self._this.reset(new CyComponent())

    property component_id:
        def __get__(self):
            return deref(self._this).component_id

    property n_merges:
        def __get__(self):
            return deref(self._this).get_n_merges()

    def __len__(self):
        return deref(self._this).get_n_tags()

    def __iter__(self):
        it = deref(self._this).tags.begin()
        while it != deref(self._this).tags.end():
            yield deref(it)
            inc(it)


cdef Component build_component(ComponentPtr ptr):
    cdef Component comp = Component()
    comp._this.reset(ptr.get())
    return comp


cdef class StreamingPartitioner:

    cdef unique_ptr[CyStreamingPartitioner] _this
    cdef weak_ptr[ComponentPtrSet] _components
    cdef weak_ptr[GuardedKmerCompMap] _tag_component_map

    cdef CyHashtable * _graph_ptr

    def __cinit__(self, graph):
        if not (isinstance(graph, Countgraph) or isinstance(graph, Nodegraph)):
            raise ValueError('Must take an object with Hashtable *')
        
        
        cdef CyCpHashtable_Object* ptr = <CyCpHashtable_Object*> graph
        self._graph_ptr = deref(ptr).hashtable
        
        self._this.reset(new CyStreamingPartitioner(self._graph_ptr))

    def consume_sequence(self, sequence):
        deref(self._this).consume_sequence(sequence.encode('utf-8'))

    def get_tag_component(self, kmer):
        cdef ComponentPtr comp
        comp = deref(self._this).get_tag_component(kmer.encode('utf-8'))
        if comp == NULL:
            return None
        else:
            return build_component(comp)

    def get_nearest_component(self, kmer):
        cdef ComponentPtr comp
        comp = deref(self._this).get_nearest_component(kmer.encode('utf-8'))
        if comp == NULL:
            return None
        else:
            return build_component(comp)

    property n_components:
        def __get__(self):
            return deref(self._this).get_n_components()
        
