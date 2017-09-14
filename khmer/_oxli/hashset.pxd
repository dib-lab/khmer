# -*- coding: UTF-8 -*-

from libcpp.set cimport set
from cython.operator cimport dereference as deref

from khmer._oxli.utils cimport _bstring
from khmer._oxli.oxli_types cimport *


cdef class HashSet:
    cdef set[HashIntoType] hs
    cdef public WordLength ksize

    cpdef add(self, HashIntoType h)
    cpdef remove(self, HashIntoType h)
