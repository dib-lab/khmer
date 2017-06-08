# cython: c_string_type=unicode, c_string_encoding=utf8
from cython.operator cimport dereference as deref

from libcpp.memory cimport unique_ptr

from .._khmer import Countgraph, Nodegraph, GraphLabels


cdef CpHashgraph * get_hashgraph_ptr(object graph):
    if not (isinstance(graph, Countgraph) or isinstance(graph, Nodegraph)):
        return NULL

    cdef CPyHashgraph_Object* ptr = <CPyHashgraph_Object*> graph
    return deref(ptr).hashgraph


cdef CpLabelHash * get_labelhash_ptr(object labels):
    if not isinstance(labels, GraphLabels):
        return NULL

    cdef CPyGraphLabels_Object * ptr = <CPyGraphLabels_Object*> labels
    return deref(ptr).labelhash


cdef class QFCounttable_t:
    cdef unique_ptr[CpQFCounttable] c_table
    def __cinit__(self, int k, int size):
        self.c_table.reset(new CpQFCounttable(k, size))

    def add(self, kmer):
        if not isinstance(kmer, unicode):
            raise ValueError("requires text input, got %s" % type(kmer))
        bytes = kmer.encode("ASCII")
        cdef const char* data = bytes
        return deref(self.c_table).add(data)

    def get_count(self, kmer):
        if not isinstance(kmer, unicode):
            raise ValueError("requires text input, got %s" % type(kmer))
        bytes = kmer.encode("ASCII")
        cdef const char* data = bytes
        return deref(self.c_table).get_count(data)
