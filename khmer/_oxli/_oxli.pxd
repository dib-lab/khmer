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
    ctypedef bool (*KmerFilter) (Kmer kmer)


cdef extern from "kmer_hash.hh" namespace "khmer":
    cdef cppclass Kmer "khmer::Kmer":
        HashIntoType kmer_f
        HashIntoType kmer_r
        HashIntoType kmer_u

        Kmer(HashIntoType, HashIntoType, HashIntoType)
        Kmer(string, WordLength)
        Kmer()

        bool is_forward() const
        void set_from_unique_hash(HashIntoType, WordLength)

    HashIntoType _hash(const string, const WordLength)
    string _revhash(HashIntoType, WordLength)


# Definitions for interacting with the CPython interface.
# All we really need are the PyObject struct definitions
# for our extension objects.
cdef extern from "_khmer.hh":
    ctypedef struct CpHashtable_Object "khmer::khmer_KHashtable_Object":
        Hashtable* hashtable


cdef extern from "hashtable.hh" namespace "khmer":
    cdef cppclass Hashtable "khmer::Hashtable":
        Hashtable(WordLength)
        const BoundedCounterType get_count(const char *) const
        const BoundedCounterType get_count(HashIntoType) const


cdef extern from "traversal.hh":
    cdef cppclass Traverser "khmer::Traverser":
        Traverser(Hashtable *)

        void push_filter(KmerFilter)
        KmerFilter pop_filter()
    
        uint32_t traverse(const Kmer&, KmerQueue&) const
        uint32_t traverse_left(const Kmer&, KmerQueue&) const     
        uint32_t traverse_right(const Kmer&, KmerQueue&) const

        uint32_t degree(const Kmer&) const
        uint32_t degree_left(const Kmer&) const
        uint32_t degree_right(const Kmer&) const





cdef extern from "partitioning.hh" namespace "khmer":

    cdef cppclass Component "khmer::Component":
        const uint64_t component_id
        set[HashIntoType] tags

        void merge(set[shared_ptr[Component]])
        void add_tag(HashIntoType)
        void add_tag(set[HashIntoType])
        uint64_t get_n_tags() const
        uint64_t get_n_created() const
        uint64_t get_n_destroyed() const

    ctypedef shared_ptr[Component] ComponentPtr
    ctypedef set[ComponentPtr] ComponentPtrSet

    cdef cppclass GuardedKmerCompMap "khmer::GuardedKmerCompMap":
        map[HashIntoType, ComponentPtr] data

        ComponentPtr get(HashIntoType)
        void set(HashIntoType, ComponentPtr)
        bool contains(HashIntoType)

    cdef cppclass StreamingPartitioner "khmer::StreamingPartitioner":
        StreamingPartitioner(Hashtable * ) except +MemoryError

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

####

