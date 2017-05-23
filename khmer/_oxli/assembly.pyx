# -*- coding: UTF-8 -*-
# cython: c_string_type=unicode, c_string_encoding=utf8

from cython.operator cimport dereference as deref

from hashing cimport CpKmer
from wrapper cimport CpLinearAssembler


cdef class LinearAssembler:

    def __cinit__(self, graph, stop_filter=None):
        self._graph_ptr = get_hashgraph_ptr(graph)
        if self._graph_ptr == NULL:
            raise ValueError('Must take an object with Hashgraph *')
        self.graph = graph
        self.set_stop_filter(stop_filter=stop_filter)
        
        if type(self) is LinearAssembler:
            self._this.reset(new CpLinearAssembler(self._graph_ptr))

    def set_stop_filter(self, stop_filter=None):
        self.stop_filter = stop_filter
        if stop_filter is not None:
            self._stop_filter_ptr = get_hashgraph_ptr(stop_filter)
            if self._stop_filter_ptr == NULL:
                raise ValueError('Must take an object with Hashgraph *')
        else:
            self._stop_filter_ptr = NULL

    cdef str _assemble(self, Kmer kmer):
        if self.stop_filter is None:
            return deref(self._this).assemble(deref(kmer._this))
        else:
            return deref(self._this).assemble(deref(kmer._this), self._stop_filter_ptr)

    def assemble(self, seed):
        if isinstance(seed, Kmer):
            return self._assemble(seed)
        else:
            return self._assemble(Kmer(str(seed)))


    cdef str _assemble_left(self, Kmer kmer):
        if self.stop_filter is None:
            return deref(self._this).assemble_left(deref(kmer._this))
        else:
            return deref(self._this).assemble_left(deref(kmer._this), self._stop_filter_ptr)

    def assemble_left(self, seed):
        if isinstance(seed, Kmer):
            return self._assemble_left(seed)
        else:
            return self._assemble_left(Kmer(str(seed)))


    cdef str _assemble_right(self, Kmer kmer):
        if self.stop_filter is None:
            return deref(self._this).assemble_right(deref(kmer._this))
        else:
            return deref(self._this).assemble_right(deref(kmer._this), self._stop_filter_ptr)

    def assemble_right(self, seed):
        if isinstance(seed, Kmer):
            return self._assemble_right(seed)
        else:
            return self._assemble_right(Kmer(str(seed)))


