from khmer._oxli.partitioning cimport StreamingPartitioner
from khmer._oxli.graphs cimport Hashgraph

cdef class PartitioningApp:

    cdef object args
    cdef readonly Hashgraph graph
    cdef readonly StreamingPartitioner partitioner
