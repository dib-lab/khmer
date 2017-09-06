# -*- coding: UTF-8 -*-
from libcpp.vector cimport vector
from libc.stdint cimport uint32_t, uint64_t
from libcpp cimport bool


cdef extern from "oxli_exception_convert.hh":
    cdef void oxli_raise_py_error()


cdef extern from "oxli/hashtable.hh" namespace "oxli":
    cdef bool _is_prime "oxli::is_prime" (uint64_t n)
    cdef vector[uint64_t] _get_n_primes_near_x "oxli::get_n_primes_near_x" (uint32_t, uint64_t)

cdef bytes _bstring(s)

cdef unicode _ustring(s)

cpdef bool is_str(object s)
cpdef bool is_num(object n)

cdef void _flatten_fill(double * fill_to, object fill_from)
cdef void _fill(double * fill_to, object fill_from)
