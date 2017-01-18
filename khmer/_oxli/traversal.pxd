from libcpp.memory cimport unique_ptr

from _oxli cimport CpTraverser, CpHashtable

cdef class Traverser:
    cdef unique_ptr[CpTraverser] _this
    cdef CpHashtable * _graph_ptr



