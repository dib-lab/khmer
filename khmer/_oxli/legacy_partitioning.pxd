from libcpp cimport bool
from libcpp.string cimport string
from libcpp.set cimport set
from libcpp.memory cimport shared_ptr, weak_ptr
from libc.stdint cimport uintptr_t as size_t

from khmer._oxli.oxli_types cimport *
from khmer._oxli.hashing cimport CpKmer
from khmer._oxli.graphs cimport CpHashgraph, CpCountgraph, Hashgraph, Countgraph
from khmer._oxli.utils cimport oxli_raise_py_error


cdef extern from "oxli/subset.hh" nogil:

    cdef struct cp_pre_partition_info "oxli::pre_partition_info":
        HashIntoType kmer
        set[HashIntoType] tagged_kmers
        cp_pre_partition_info(HashIntoType)

    cdef cppclass CpSubsetPartition "oxli::SubsetPartition":
        CpSubsetPartition(CpHashgraph *)
        
        PartitionID assign_partition_id(HashIntoType, HashIntoTypeSet &)
        void set_partition_id(HashIntoType, PartitionID)
        void set_partition_id(string, PartitionID)
        PartitionID join_partitions(PartitionID, PartitionID)
        PartitionID get_partition_id(string)
        PartitionID get_partition_id(HashIntoType)

        PartitionID * get_new_partition()
        void merge(CpSubsetPartition *)
        void merge_from_disk(string) except +oxli_raise_py_error
        void _merge_from_disk_consolidate(PartitionPtrMap&)

        void save_partitionmap(string) except +oxli_raise_py_error 
        void load_partitionmap(string) except +oxli_raise_py_error
        void _validate_pmap()
        
        void find_all_tags(CpKmer, HashIntoTypeSet &, HashIntoTypeSet &,
                           bool, bool) except +oxli_raise_py_error
        void find_all_tags(CpKmer, HashIntoTypeSet &, 
                           HashIntoTypeSet &) except +oxli_raise_py_error

        unsigned int sweep_for_tags(const string&, HashIntoTypeSet &,
                                    HashIntoTypeSet &, unsigned int,
                                    bool, bool)

        void find_all_tags_truncate_on_abundance(CpKmer, HashIntoTypeSet &,
                                                 HashIntoTypeSet &,
                                                 BoundedCounterType,
                                                 BoundedCounterType,
                                                 bool,
                                                 bool) except +oxli_raise_py_error

        void do_partition(HashIntoType, HashIntoType, 
                          bool, bool) except +oxli_raise_py_error

        void do_partition_with_abundance(HashIntoType, HashIntoType,
                                         BoundedCounterType,
                                         BoundedCounterType,
                                         bool, bool) except +oxli_raise_py_error

        void count_partitions(size_t&, size_t&)

        size_t output_partitioned_file(string &, 
                                       string &, 
                                       bool) except +oxli_raise_py_error

        void partition_sizes(PartitionCountMap&, unsigned int &) const
        void partition_size_distribution(PartitionCountDistribution &,
                                         unsigned int &) const
        void partition_average_coverages(PartitionCountMap &,
                                         CpCountgraph *) const
        unsigned long long repartition_largest_partition(unsigned int,
                                                         unsigned int,
                                                         unsigned int,
                                                         CpCountgraph&)
        void repartition_a_partition(const HashIntoTypeSet &) except +oxli_raise_py_error
        void _clear_partition(PartitionID, HashIntoTypeSet &)
        void _merge_other(HashIntoType, PartitionID, PartitionPtrMap &)
        void report_on_partitions()
        
cdef class PrePartitionInfo:
    cdef shared_ptr[cp_pre_partition_info] _this
    @staticmethod
    cdef PrePartitionInfo wrap(cp_pre_partition_info *)
    @staticmethod
    cdef PrePartitionInfo create(CpKmer kmer)


cdef class SubsetPartition:
    cdef shared_ptr[CpSubsetPartition] _this
