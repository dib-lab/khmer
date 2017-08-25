# -*- coding: UTF-8 -*-
# cython: c_string_type=unicode, c_string_encoding=utf8

from cython.operator cimport dereference as deref

from utils cimport _bstring


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

    cdef str _assemble(self, Kmer kmer):
        if self.stop_filter is None:
            return deref(self._this).assemble(deref(kmer._this))
        else:
            return deref(self._this).assemble(deref(kmer._this), 
                                              self._stop_filter_ptr.get())

    def assemble(self, seed):
        if isinstance(seed, Kmer):
            return self._assemble(seed)
        else:
            return self._assemble(Kmer(str(seed)))


    cdef str _assemble_left(self, Kmer kmer):
        if self.stop_filter is None:
            return deref(self._this).assemble_left(deref(kmer._this))
        else:
            return deref(self._this).assemble_left(deref(kmer._this), 
                                                   self._stop_filter_ptr.get())

    def assemble_left(self, seed):
        if isinstance(seed, Kmer):
            return self._assemble_left(seed)
        else:
            return self._assemble_left(Kmer(str(seed)))


    cdef str _assemble_right(self, Kmer kmer):
        if self.stop_filter is None:
            return deref(self._this).assemble_right(deref(kmer._this))
        else:
            return deref(self._this).assemble_right(deref(kmer._this), 
                                                    self._stop_filter_ptr.get())

    def assemble_right(self, seed):
        if isinstance(seed, Kmer):
            return self._assemble_right(seed)
        else:
            return self._assemble_right(Kmer(str(seed)))


cdef class SimpleLabeledAssembler:

    def __cinit__(self, GraphLabels labels not None, Hashgraph stop_filter=None):
        self.labels = labels
        self._label_ptr = labels._lh_this
        self.set_stop_filter(stop_filter=stop_filter)
        
        if type(self) is SimpleLabeledAssembler:
            self._this.reset(new CpSimpleLabeledAssembler(self._label_ptr.get()))

    def set_stop_filter(self, Hashgraph stop_filter=None):
        self.stop_filter = stop_filter
        if stop_filter is not None:
            self._stop_filter_ptr = stop_filter._hg_this

    cdef vector[string] _assemble(self, Kmer kmer):
        if self.stop_filter is None:
            return deref(self._this).assemble(deref(kmer._this))
        else:
            return deref(self._this).assemble(deref(kmer._this),
                                              self._stop_filter_ptr.get())

    def assemble(self, seed):
        if isinstance(seed, Kmer):
            return self._assemble(seed)
        else:
            return self._assemble(Kmer(str(seed)))


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

    cdef vector[string] _assemble(self, Kmer kmer):
        if self.stop_filter is None:
            return deref(self._this).assemble(deref(kmer._this))
        else:
            return deref(self._this).assemble(deref(kmer._this),
                                              self._stop_filter_ptr.get())

    def consume(self, sequence):
        return deref(self._this).consume(_bstring(sequence))

    def assemble(self, seed):
        if isinstance(seed, Kmer):
            return self._assemble(seed)
        else:
            return self._assemble(Kmer(str(seed)))
