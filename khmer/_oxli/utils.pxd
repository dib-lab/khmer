# -*- coding: UTF-8 -*-
from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdint cimport uint32_t, uint64_t
from libcpp cimport bool


cdef extern from "oxli/oxli_exception_convert.hh":
    cdef void oxli_raise_py_error()


cdef extern from "oxli/hashtable.hh" namespace "oxli":
    cdef bool _is_prime "oxli::is_prime" (uint64_t n)
    cdef vector[uint64_t] _get_n_primes_near_x "oxli::get_n_primes_near_x" (uint32_t, uint64_t)

cdef extern from "oxli/oxli.hh":
    cdef string _get_version_cpp "oxli::get_version_cpp" ()
    cdef const char * SAVED_SIGNATURE
    cdef int SAVED_FORMAT_VERSION
    cdef int SAVED_COUNTING_HT
    cdef int SAVED_HASHBITS
    cdef int SAVED_TAGS
    cdef int SAVED_STOPTAGS
    cdef int SAVED_SUBSET
    cdef int SAVED_LABELSET
    cdef int SAVED_SMALLCOUNT
    cdef int SAVED_QFCOUNT


cdef bytes _bstring(s)

cdef unicode _ustring(s)

cpdef bool is_str(object s)
cpdef bool is_num(object n)

cdef void _flatten_fill(double * fill_to, object fill_from)
cdef void _fill(double * fill_to, object fill_from)

cpdef str get_version_cpp()
