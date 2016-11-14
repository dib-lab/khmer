from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.queue cimport queue
from libcpp.memory cimport unique_ptr, weak_ptr, shared_ptr
from libcpp.utility cimport pair
from libc.stdint cimport uint32_t, uint8_t, uint64_t


cdef extern from "khmer.hh" namespace "khmer":
    ctypedef unsigned long long int HashIntoType
    ctypedef unsigned char WordLength
    ctypedef unsigned short int BoundedCounterType
    ctypedef queue[Kmer] KmerQueue


cdef extern from "hashtable.hh" namespace "khmer":
    cdef cppclass CyHashtable "khmer::Hashtable":
        CyHashtable(WordLength)
        const BoundedCounterType get_count(const char *) const
        const BoundedCounterType get_count(HashIntoType) const


cdef extern from "_khmer.hh":
    ctypedef struct CyCpHashtable_Object "khmer::khmer_KHashtable_Object":
        CyHashtable* hashtable


cdef extern from "kmer_hash.hh" namespace "khmer":
    cdef cppclass Kmer:
        HashIntoType kmer_f
        HashIntoType kmer_r
        HashIntoType kmer_u

        Kmer(HashIntoType, HashIntoType, HashIntoType)
        Kmer(string, WordLength)
        Kmer()

        bool is_forward() const




cdef extern from "traversal.hh":
    cdef cppclass CyLeftNodeGatherer "khmer::NodeGatherer<0>":
        CyLeftNodeGatherer(CyHashtable *)
        void push_filter()


cdef extern from "partitioning.hh" namespace "khmer":

    cdef cppclass CyComponent "khmer::Component":
        const uint64_t component_id
        set[HashIntoType] tags

        void merge(set[shared_ptr[CyComponent]])
        void add_tag(HashIntoType)
        void add_tag(set[HashIntoType])
        uint64_t get_n_tags() const
        uint64_t get_n_created() const
        uint64_t get_n_destroyed() const

    ctypedef shared_ptr[CyComponent] ComponentPtr
    ctypedef set[ComponentPtr] ComponentPtrSet

    cdef cppclass GuardedKmerCompMap "khmer::GuardedKmerCompMap":
        map[HashIntoType, ComponentPtr] data

        ComponentPtr get(HashIntoType)
        void set(HashIntoType, ComponentPtr)
        bool contains(HashIntoType)

    cdef cppclass CyStreamingPartitioner "khmer::StreamingPartitioner":
        CyStreamingPartitioner(CyHashtable * ) except +MemoryError

        void consume_sequence(string&) except +MemoryError
        uint64_t consume_fasta(string&) except +MemoryError
        void map_tags_to_component(set[HashIntoType]&, ComponentPtr&)
        void find_connected_tags(queue[Kmer]&, 
                                 set[HashIntoType]&,
                                 set[HashIntoType]&) except +MemoryError
        uint64_t get_n_components() const
        uint64_t get_n_tags() const

        ComponentPtr get_tag_component(string&) const
        ComponentPtr get_nearest_component(string&) const

        weak_ptr[ComponentPtrSet] get_component_set()
        weak_ptr[GuardedKmerCompMap] get_tag_component_map()
