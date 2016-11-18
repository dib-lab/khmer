from partitioning cimport StreamingPartitioner

cdef class PartitioningApp:

    cdef object args
    cdef object graph
    cdef StreamingPartitioner partitioner
