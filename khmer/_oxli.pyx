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


cdef class StreamingPartitioner:

    cdef unique_ptr[CyStreamingPartitioner] _this
    cdef CyHashtable * _graph_ptr

    def __cinit__(self, graph):
        if not (isinstance(graph, Countgraph) or isinstance(graph, Nodegraph)):
            raise ValueError('Must take an object with Hashtable *')
        
        
        cdef CyCpHashtable_Object* ptr = <CyCpHashtable_Object*> graph
        self._graph_ptr = deref(ptr).hashtable
        
        self._this.reset(new CyStreamingPartitioner(self._graph_ptr))

    def consume_sequence(self, sequence):
        deref(self._this).consume_sequence(sequence.encode('utf-8'))

    property n_components:
        def __get__(self):
            return deref(self._this).get_n_components()
        
