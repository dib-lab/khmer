from libcpp.memory cimport unique_ptr
from cython.operator cimport dereference as deref

from .._khmer import Countgraph
from .._khmer import Nodegraph

from _oxli cimport *
cimport hashing
import hashing

cdef class Traverser:

    def __cinit__(self, graph):
        if not (isinstance(graph, Countgraph) or isinstance(graph, Nodegraph)):
            raise ValueError('Must take an object with Hashtable *')
        
        cdef CPyHashtable_Object* ptr = <CPyHashtable_Object*> graph
        self._graph_ptr = deref(ptr).hashtable
        
        self._this.reset(new CpTraverser(self._graph_ptr))
    '''
    def left_neighbors(self, Kmer kmer):
        cdef KmerQueue kmer_q
        cdef Kmer neighbor
        deref(self._this).traverse_left(kmer, kmer_q)
        for neighbor in kmer_q:
            yield neighbor
    '''

    def right_neighbors(self, hashing.Kmer kmer):
        cdef KmerQueue kmer_q
        cdef CpKmer neighbor
        cdef CpKmer* ptr
        deref(self._this).traverse_right(deref(kmer._this), kmer_q)
        while(kmer_q.empty() == 0):
            neighbor = kmer_q.back()
            ptr = new CpKmer(neighbor)
            yield hashing.Kmer.wrap(ptr, deref(self._graph_ptr).ksize())
            kmer_q.pop()
