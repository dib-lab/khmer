import cython
from cython.operator cimport dereference as deref, preincrement as inc
from libc.limits cimport UINT_MAX
from libcpp.memory cimport unique_ptr, weak_ptr
from libc.stdint cimport uintptr_t
from libc.stdio cimport FILE, fopen, fwrite, fclose, stdout, stderr, fprintf

from khmer._khmer import Countgraph
from khmer._khmer import Nodegraph

cdef class Component:

    cdef ComponentPtr _this

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
        it = deref(self._this).tags.begin()
        while it != deref(self._this).tags.end():
            yield deref(it)
            inc(it)

    def __hash__(self):
        return <uintptr_t>self._this.get()

    def __richcmp__(x, y, op):
        if op == 2:
            return x.component_id == y.component_id
        else:
            raise NotImplementedError('Operator not available.')

    @staticmethod
    cdef Component create(ComponentPtr ptr):
        cdef Component comp = Component()
        comp._this = ptr
        return comp

    @staticmethod
    cdef vector[BoundedCounterType] _tag_counts(ComponentPtr comp, CyHashtable * graph):
        cdef uint64_t n_tags = deref(comp).get_n_tags()
        cdef vector[BoundedCounterType] counts
        counts = vector[BoundedCounterType](n_tags)
        cdef int idx
        cdef uint64_t tag
        for idx, tag in deref(comp).tags:
            counts[idx] = deref(graph).get_count(tag)
        return counts

    @staticmethod
    cdef float _mean_tag_count(ComponentPtr comp, CyHashtable * graph):
        cdef uint64_t n_tags = deref(comp).get_n_tags()
        cdef float acc = 0
        cdef uint64_t tag
        for tag in deref(comp).tags:
            acc += <float>deref(graph).get_count(tag)
        return acc / <float>n_tags


cdef class StreamingPartitioner:

    cdef unique_ptr[CyStreamingPartitioner] _this
    cdef weak_ptr[ComponentPtrSet] _components
    cdef weak_ptr[GuardedKmerCompMap] _tag_component_map
    cdef CyHashtable * _graph_ptr
    cdef readonly uint64_t n_consumed

    def __cinit__(self, graph):
        if not (isinstance(graph, Countgraph) or isinstance(graph, Nodegraph)):
            raise ValueError('Must take an object with Hashtable *')
        
        
        cdef CyCpHashtable_Object* ptr = <CyCpHashtable_Object*> graph
        self._graph_ptr = deref(ptr).hashtable
        
        self._this.reset(new CyStreamingPartitioner(self._graph_ptr))

        self._tag_component_map = deref(self._this).get_tag_component_map()
        self._components = deref(self._this).get_component_set()
        self.n_consumed = 0

    def consume_sequence(self, sequence):
        deref(self._this).consume_sequence(sequence.encode('utf-8'))
        self.n_consumed += 1

    def get_tag_component(self, kmer):
        cdef ComponentPtr compptr
        compptr = deref(self._this).get_tag_component(kmer.encode('utf-8'))
        if compptr == NULL:
            return None
        else:
            return Component.create(compptr)

    def get_nearest_component(self, kmer):
        cdef ComponentPtr compptr
        compptr = deref(self._this).get_nearest_component(kmer.encode('utf-8'))
        if compptr == NULL:
            return None
        else:
            return Component.create(compptr)

    def components(self):
        cdef shared_ptr[ComponentPtrSet] locked
        lockedptr = self._components.lock()
        if lockedptr:
            it = deref(lockedptr).begin()
            while it != deref(lockedptr).end():
                yield Component.create(deref(it))
                inc(it)
        else:
            raise MemoryError("Can't locked underlying Component set")

    def tag_components(self):
        cdef shared_ptr[GuardedKmerCompMap] locked
        lockedptr = self._tag_component_map.lock()
        if lockedptr:
            it = deref(lockedptr).data.begin()
            while it != deref(lockedptr).data.end():
                yield deref(it).first, Component.create(deref(it).second)
                inc(it)
        else:
            raise MemoryError("Can't locked underlying Component set")

    def write_components(self, filename):
        cdef FILE* fp
        fp = fopen(filename.encode('utf-8'), 'wb')
        if fp == NULL:
            raise IOError("Can't open file.")
        
        cdef ComponentPtr cmpptr
        cdef shared_ptr[ComponentPtrSet] locked
        lockedptr = self._components.lock()

        if lockedptr:      
            it = deref(lockedptr).begin()
            while it != deref(lockedptr).end():
                cmpptr = deref(it)

                fprintf(fp, "%u,%u,%f\n", 
                        deref(cmpptr).component_id,
                        deref(cmpptr).get_n_tags(),
                        Component._mean_tag_count(cmpptr, self._graph_ptr))
                inc(it)
        fclose(fp)

    property n_components:
        def __get__(self):
            return deref(self._this).get_n_components()
        
    property n_tags:
        def __get__(self):
            return deref(self._this).get_n_tags()
