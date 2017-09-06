# -*- coding: UTF-8 -*-

from cython.operator cimport dereference as deref

from khmer._oxli.utils cimport _bstring


cdef class LinearAssembler:

    def __cinit__(self, Hashgraph graph not None, Hashgraph stop_filter=None):
        self.graph = graph
        self._graph_ptr = graph._hg_this
        self.set_stop_filter(stop_filter=stop_filter)
        
        if type(self) is LinearAssembler:
            self._this = make_shared[CpLinearAssembler](self._graph_ptr.get())

    def set_stop_filter(self, Hashgraph stop_filter=None):
        self.stop_filter = stop_filter
        if stop_filter is not None:
            self._stop_filter_ptr = stop_filter._hg_this

    cdef str _assemble(self, CpKmer kmer):
        if self.stop_filter is None:
            return deref(self._this).assemble(kmer)
        else:
            return deref(self._this).assemble(kmer, 
                                              self._stop_filter_ptr.get())

    def assemble(self, object seed):
        cdef CpKmer _seed = self.graph._build_kmer(seed)
        return self._assemble(_seed)

    cdef str _assemble_left(self, CpKmer kmer):
        if self.stop_filter is None:
            return deref(self._this).assemble_left(kmer)
        else:
            return deref(self._this).assemble_left(kmer, 
                                                   self._stop_filter_ptr.get())

    def assemble_left(self, object seed):
        cdef CpKmer _seed = self.graph._build_kmer(seed)
        return self._assemble_left(_seed)

    cdef str _assemble_right(self, CpKmer kmer):
        if self.stop_filter is None:
            return deref(self._this).assemble_right(kmer)
        else:
            return deref(self._this).assemble_right(kmer, 
                                                    self._stop_filter_ptr.get())

    def assemble_right(self, object seed):
        cdef CpKmer _seed = self.graph._build_kmer(seed)
        return self._assemble_right(_seed)


cdef class SimpleLabeledAssembler:

    def __cinit__(self, GraphLabels labels not None, Hashgraph stop_filter=None):
        self.labels = labels
        self.graph = labels.graph
        self._label_ptr = labels._lh_this
        self.set_stop_filter(stop_filter=stop_filter)
        
        if type(self) is SimpleLabeledAssembler:
            self._this.reset(new CpSimpleLabeledAssembler(self._label_ptr.get()))

    def set_stop_filter(self, Hashgraph stop_filter=None):
        self.stop_filter = stop_filter
        if stop_filter is not None:
            self._stop_filter_ptr = stop_filter._hg_this

    cdef vector[string] _assemble(self, CpKmer kmer):
        if self.stop_filter is None:
            return deref(self._this).assemble(kmer)
        else:
            return deref(self._this).assemble(kmer,
                                              self._stop_filter_ptr.get())

    def assemble(self, object seed):
        cdef CpKmer _seed = self.graph._build_kmer(seed)
        return self._assemble(_seed)


cdef class JunctionCountAssembler:

    def __cinit__(self, Hashgraph graph not None, Hashgraph stop_filter=None):
        self.graph = graph
        self._graph_ptr = graph._hg_this
        self.set_stop_filter(stop_filter=stop_filter)
        
        if type(self) is JunctionCountAssembler:
            self._this = make_shared[CpJunctionCountAssembler](self._graph_ptr.get())

    def set_stop_filter(self, Hashgraph stop_filter=None):
        self.stop_filter = stop_filter
        if stop_filter is not None:
            self._stop_filter_ptr = stop_filter._hg_this

    cdef vector[string] _assemble(self, CpKmer kmer):
        if self.stop_filter is None:
            return deref(self._this).assemble(kmer)
        else:
            return deref(self._this).assemble(kmer,
                                              self._stop_filter_ptr.get())

    def consume(self, str sequence):
        return deref(self._this).consume(_bstring(sequence))

    def assemble(self, object seed):
        cdef CpKmer _seed = self.graph._build_kmer(seed)
        return self._assemble(_seed)
