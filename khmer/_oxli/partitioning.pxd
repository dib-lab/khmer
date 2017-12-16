from libcpp cimport bool
from libcpp.memory cimport unique_ptr, weak_ptr, shared_ptr
from libcpp.unordered_map cimport unordered_map
from libcpp.vector cimport vector
from libcpp.set cimport set
from libcpp.queue cimport queue
from libcpp.string cimport string
from libc.stdint cimport uint32_t, uint8_t, uint64_t
from libc.stdio cimport FILE

from khmer._oxli.hashing cimport CpKmer, Kmer, KmerQueue
from khmer._oxli.hist cimport CpHistogram
from khmer._oxli.graphs cimport CpHashgraph, Hashgraph
from khmer._oxli.oxli_types cimport *


cdef extern from "oxli/partitioning.hh" namespace "oxli":

    ctypedef vector[HashIntoType] TagVector

    cdef cppclass CpComponent "oxli::Component":
        CpComponent()
        CpComponent(uint64_t)

        CpHistogram coverage
        const uint64_t component_id
        vector[HashIntoType] tags

        void kill()
        bool is_alive() const

        void add_tag(HashIntoType)
        void add_tags(TagVector&)

        uint64_t get_n_tags() const
        uint64_t get_n_created() const
        uint64_t get_n_destroyed() const

        void update_coverage(CpHashgraph *) 

    ctypedef shared_ptr[CpComponent] ComponentPtr
    ctypedef set[ComponentPtr] ComponentPtrSet
    ctypedef vector[ComponentPtr] ComponentPtrVector

    cdef cppclass CpGuardedHashCompMap "oxli::GuardedHashCompMap":
        unordered_map[HashIntoType, ComponentPtr] data

        ComponentPtr get(HashIntoType)
        void set(HashIntoType, ComponentPtr)
        bool contains(HashIntoType)

    cdef cppclass CpComponentMap "oxli::ComponentMap":
        CpComponentMap(WordLength, WordLength, uint64_t)

        void create_component(TagVector&)
        uint32_t create_and_merge_components(TagVector&)
        void map_tags_to_component(TagVector&, ComponentPtr&)
        uint32_t merge_components(ComponentPtr&, ComponentPtrSet&)

        bool contains(HashIntoType)
        ComponentPtr get(HashIntoType) const

        uint64_t get_n_components() const
        uint64_t get_n_tags() const
        weak_ptr[ComponentPtrVector] get_components()
        weak_ptr[CpGuardedHashCompMap] get_tag_component_map()

    cdef cppclass CpStreamingPartitioner "oxli::StreamingPartitioner" (CpComponentMap):
        CpStreamingPartitioner(CpHashgraph * ) except +MemoryError
        CpStreamingPartitioner(CpHashgraph *, uint32_t) except +MemoryError
        
        CpHashgraph * graph
        uint64_t consume(string&) nogil except +MemoryError
        uint64_t  consume_pair(string&, string&) nogil except +MemoryError
        uint64_t consume_fasta(string&) except +MemoryError

        uint64_t seed_sequence(string&, TagVector&, KmerQueue&,
                           set[HashIntoType]&) except +MemoryError

        void find_connected_tags(KmerQueue&, 
                                 TagVector&,
                                 set[HashIntoType]&) except +MemoryError

        void find_connected_tags(KmerQueue&, 
                                 TagVector&,
                                 set[HashIntoType]&,
                                 bool) except +MemoryError

        ComponentPtr find_nearest_component(string&) const
        ComponentPtr find_nearest_component(CpKmer) const

        uint64_t get_n_consumed() const
        uint32_t get_tag_density() const

        ComponentPtr get(string&) const


cdef class Component:
    cdef ComponentPtr _this

    cdef void save(self, FILE * fp)

    @staticmethod
    cdef Component wrap(ComponentPtr ptr)

    @staticmethod
    cdef vector[BoundedCounterType] _tag_counts(ComponentPtr comp, CpHashgraph* graph)

    @staticmethod
    cdef float _mean_tag_count(ComponentPtr comp, CpHashgraph * graph)

    @staticmethod
    cdef ComponentPtr load(uint64_t component_id, list tags)


cdef class StreamingPartitioner:
    cdef shared_ptr[CpStreamingPartitioner] _this
    cdef weak_ptr[ComponentPtrVector] _components
    cdef weak_ptr[CpGuardedHashCompMap] _tag_component_map
    cdef public Hashgraph graph
    cdef readonly uint64_t n_consumed

