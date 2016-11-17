import cython
from cython.operator cimport dereference as deref, preincrement as inc

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.queue cimport queue
from libcpp.memory cimport unique_ptr, weak_ptr, shared_ptr
from libcpp.utility cimport pair

from libc.stdint cimport uint32_t, uint8_t, uint64_t
from libc.limits cimport UINT_MAX

from libc.stdint cimport uintptr_t
from libc.stdio cimport FILE, fopen, fwrite, fclose, stdout, stderr, fprintf

from _oxli cimport *
from .._khmer import Countgraph
from .._khmer import Nodegraph

cdef class Component:

    def __cinit__(self, Component other=None):
        if other is not None:
            self._this.reset(other._this.get())

    property component_id:
        def __get__(self):
            return deref(self._this).component_id

    property _n_created:
        def __get__(self):
            return deref(self._this).get_n_created()

    property _n_destroyed:
        def __get__(self):
            return deref(self._this).get_n_destroyed()

    def __len__(self):
        return deref(self._this).get_n_tags()

    def __iter__(self):
        cdef HashIntoType tag
        for tag in deref(self._this).tags:
            yield tag

    def __hash__(self):
        return <uintptr_t>self._this.get()

    def __richcmp__(x, y, op):
        if op == 2:
            return x.component_id == y.component_id
        else:
            raise NotImplementedError('Operator not available.')

    @staticmethod
    cdef Component wrap(ComponentPtr ptr):
        cdef Component comp = Component()
        comp._this = ptr
        return comp

    @staticmethod
    cdef vector[BoundedCounterType] _tag_counts(ComponentPtr comp, CpHashtable * graph):
        cdef uint64_t n_tags = deref(comp).get_n_tags()
        cdef vector[BoundedCounterType] counts
        counts = vector[BoundedCounterType](n_tags)
        cdef int idx
        cdef uint64_t tag
        for idx, tag in deref(comp).tags:
            counts[idx] = deref(graph).get_count(tag)
        return counts

    @staticmethod
    def tag_counts(Component component, graph):
        cdef CPyHashtable_Object* graph_ptr = <CPyHashtable_Object*> graph
        return Component._tag_counts(component._this, deref(graph_ptr).hashtable)

    @staticmethod
    cdef float _mean_tag_count(ComponentPtr comp, CpHashtable * graph):
        cdef uint64_t n_tags = deref(comp).get_n_tags()
        cdef float acc = 0
        cdef uint64_t tag
        for tag in deref(comp).tags:
            acc += <float>deref(graph).get_count(tag)
        return acc / <float>n_tags


cdef class StreamingPartitioner:

    def __cinit__(self, graph):
        if not (isinstance(graph, Countgraph) or isinstance(graph, Nodegraph)):
            raise ValueError('Must take an object with Hashtable *')
        
        cdef CPyHashtable_Object* ptr = <CPyHashtable_Object*> graph
        self._graph_ptr = deref(ptr).hashtable
        
        self._this.reset(new CpStreamingPartitioner(self._graph_ptr))

        self._tag_component_map = deref(self._this).get_tag_component_map()
        self._components = deref(self._this).get_component_set()
        self.n_consumed = 0

    def consume(self, sequence):
        deref(self._this).consume(sequence.encode('utf-8'))
        self.n_consumed += 1

    def consume_pair(self, first, second):
        deref(self._this).consume_pair(first.encode('utf-8'),
                                       second.encode('utf-8'))
        self.n_consumed += 2

    def consume_fasta(self, filename):
        return deref(self._this).consume_fasta(filename.encode('utf-8'))

    def get_tag_component(self, kmer):
        cdef ComponentPtr compptr
        compptr = deref(self._this).get_tag_component(kmer.encode('utf-8'))
        if compptr == NULL:
            return None
        else:
            return Component.wrap(compptr)

    def get_nearest_component(self, kmer):
        cdef ComponentPtr compptr
        compptr = deref(self._this).get_nearest_component(kmer.encode('utf-8'))
        if compptr == NULL:
            return None
        else:
            return Component.wrap(compptr)

    def components(self):
        cdef shared_ptr[ComponentPtrSet] locked
        cdef ComponentPtr cmpptr
        lockedptr = self._components.lock()
        if lockedptr:
            for cmpptr in deref(lockedptr):
                yield Component.wrap(cmpptr)
        else:
            raise MemoryError("Can't locked underlying Component set")

    def tag_components(self):
        cdef shared_ptr[CpGuardedKmerCompMap] locked
        cdef pair[HashIntoType,ComponentPtr] cpair
        lockedptr = self._tag_component_map.lock()
        if lockedptr:
            for cpair in deref(lockedptr).data:
                yield cpair.first, Component.wrap(cpair.second)
        else:
            raise MemoryError("Can't locked underlying Component set")

    def write_components(self, filename):
        cdef FILE* fp
        fp = fopen(filename.encode('utf-8'), 'wb')
        if fp == NULL:
            raise IOError("Can't open file.")
        
        cdef ComponentPtr cmpptr
        cdef shared_ptr[ComponentPtrSet] lockedptr
        lockedptr = self._components.lock()

        if lockedptr:      
            for cmpptr in deref(lockedptr):
                fprintf(fp, "%llu,%llu,%f\n", 
                        deref(cmpptr).component_id,
                        deref(cmpptr).get_n_tags(),
                        Component._mean_tag_count(cmpptr, self._graph_ptr))
        fclose(fp)

    property n_components:
        def __get__(self):
            return deref(self._this).get_n_components()
        
    property n_tags:
        def __get__(self):
            return deref(self._this).get_n_tags()
