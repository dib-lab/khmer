from cython.operator cimport dereference as deref

from libcpp.memory cimport make_shared

from khmer._oxli.utils cimport _bstring
from khmer._oxli.graphs cimport CpHashgraph, CpCountgraph, Hashgraph, Countgraph

cdef class PrePartitionInfo:

    @staticmethod
    cdef PrePartitionInfo wrap(cp_pre_partition_info * ptr):
        cdef PrePartitionInfo part = PrePartitionInfo()
        part._this.reset(ptr)
        return part

    @staticmethod
    cdef PrePartitionInfo create(CpKmer kmer):
        cdef PrePartitionInfo part = PrePartitionInfo()
        part._this = make_shared[cp_pre_partition_info](kmer)
        return part


cdef class SubsetPartition:

    def __cinit__(self, Hashgraph graph):
        if type(self) is SubsetPartition:
            self._this = make_shared[CpSubsetPartition](graph._hg_this.get())

    def count_partitions(self):
        cdef size_t n_partitions = 0
        cdef size_t n_unassigned = 0
        deref(self._this).count_partitions(n_partitions, n_unassigned)

        return n_partitions, n_unassigned

    def report_on_partitions(self):
        deref(self._this).report_on_partitions()

    def partition_size_distribution(self):
        cdef PartitionCountDistribution d
        cdef unsigned int n_unassigned = 0
        deref(self._this).partition_size_distribution(d, n_unassigned)

        cdef list dist = []
        for pair in d:
            dist.append((pair.first, pair.second))

        return dist, n_unassigned

    def partition_sizes(self, unsigned int min_size=0):
        cdef PartitionCountMap cm
        cdef unsigned int n_unassigned = 0
        deref(self._this).partition_sizes(cm, n_unassigned)

        cdef list sizes = []
        for pair in cm:
            if pair.second >= min_size:
                sizes.append((pair.first, pair.second))

        return sizes, n_unassigned

    def partition_average_coverages(self, Countgraph graph not None):
        cdef PartitionCountMap cm
        deref(self._this).partition_average_coverages(cm, graph._cg_this.get())

        cdef list covs = []
        for pair in cm:
            covs.append((pair.first, pair.second))

        return covs

    def save_partitionmap(self, str filename):
        cdef string _filename = _bstring(filename)
        with nogil:
            deref(self._this).save_partitionmap(_filename)

    def load_partitionmap(self, str filename):
        cdef string _filename = _bstring(filename)
        with nogil:
            deref(self._this).load_partitionmap(_filename)

    @staticmethod
    def load(str filename, Hashgraph graph):
        cdef SubsetPartition subset = SubsetPartition(graph)
        subset.load_partitionmap(filename)
        return subset

    def _validate_partitionmap(self):
        deref(self._this)._validate_pmap()
