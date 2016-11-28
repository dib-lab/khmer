from partitioning cimport StreamingPartitioner

cdef class PartitioningApp:

    cdef object args
    cdef readonly object graph
    cdef readonly StreamingPartitioner partitioner
