from libcpp.vector cimport vector
from libc.stdint cimport uint32_t, uint64_t
from libcpp cimport bool


cdef extern from "hashtable.hh" namespace "khmer":
    cdef bool _is_prime "khmer::is_prime" (uint64_t n)
    cdef vector[uint64_t] _get_n_primes_near_x "khmer::get_n_primes_near_x" (uint32_t, uint64_t)
