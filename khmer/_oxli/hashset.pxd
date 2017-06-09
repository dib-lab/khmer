# -*- coding: UTF-8 -*-
# cython: c_string_type=unicode, c_string_encoding=utf8

from libcpp.set cimport set
from cython.operator cimport dereference as deref

from utils cimport _bstring

from oxli_types cimport *


cdef class HashSet:
    cdef set[HashIntoType] hs
    cdef public WordLength ksize

    cpdef add(self, HashIntoType h)
    cpdef remove(self, HashIntoType h)
