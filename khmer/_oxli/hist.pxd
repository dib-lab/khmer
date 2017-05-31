from libc.stdint cimport uint64_t

cdef extern from "oxli/hist.hh" namespace "oxli":

    cdef cppclass CpHistogram "oxli::Histogram<16>":
        uint64_t[16] bins

        CpHistogram()

        void add(uint64_t)
        void clear()
