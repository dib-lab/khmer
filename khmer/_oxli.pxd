from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.queue cimport queue
from libcpp.memory cimport unique_ptr, weak_ptr, shared_ptr
from libc.stdint cimport uint32_t, uint8_t, uint64_t


cdef extern from "khmer.hh" namespace "khmer":
    ctypedef unsigned long long int HashIntoType
    ctypedef unsigned char WordLength
    ctypedef unsigned short int BoundedCounterType
    ctypedef queue[Kmer] KmerQueue


cdef extern from "hashtable.hh" namespace "khmer":
    cdef cppclass CyHashtable "khmer::Hashtable":
        CyHashtable(WordLength)


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


cdef extern from "partitioning.hh" namespace "khmer":
    cdef cppclass CyComponent "khmer::Component":
        const uint64_t component_id
        set[HashIntoType] tags

        void merge(set[shared_ptr[CyComponent]])
        void add_tag(HashIntoType)
        void add_tag(set[HashIntoType])
        uint64_t get_n_tags()
        uint64_t get_n_merges()

    ctypedef shared_ptr[CyComponent] ComponentPtr


cdef extern from "partitioning.hh" namespace "khmer":
    cdef cppclass CyStreamingPartitioner "khmer::StreamingPartitioner":
        CyStreamingPartitioner(CyHashtable * ) except +MemoryError

        void consume_sequence(string&) except +MemoryError
        void map_tags_to_component(set[HashIntoType], ComponentPtr)
        void find_connected_tags(queue[Kmer]&, 
                                 set[HashIntoType]&,
                                 set[HashIntoType]&) except +MemoryError
        uint64_t get_n_components()
