from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.queue cimport queue
from libcpp.memory cimport unique_ptr, weak_ptr, shared_ptr
from libcpp.utility cimport pair
from libc.stdint cimport uint32_t, uint8_t, uint64_t


########################################################################
#
# Core: typedefs from khmer.hh.
#
########################################################################

cdef extern from "khmer.hh" namespace "khmer":
    ctypedef unsigned long long int HashIntoType
    ctypedef unsigned char WordLength
    ctypedef unsigned short int BoundedCounterType
    ctypedef queue[CpKmer] KmerQueue
    ctypedef bool (*KmerFilter) (CpKmer kmer)


########################################################################
#
# Hashing: Definitions from kmer_hash.hh.
#
########################################################################

cdef extern from "kmer_hash.hh" namespace "khmer":
    cdef cppclass CpKmer "khmer::Kmer":
        HashIntoType kmer_f
        HashIntoType kmer_r
        HashIntoType kmer_u

        CpKmer(HashIntoType, HashIntoType, HashIntoType)
        CpKmer(string, WordLength)
        CpKmer(const CpKmer&)
        CpKmer()

        bool is_forward() const
        void set_from_unique_hash(HashIntoType, WordLength)

    HashIntoType _hash(const string, const WordLength)
    string _revhash(HashIntoType, WordLength)


########################################################################
#
# ReadParser: read parsing stuff, Read object
#
########################################################################


cdef extern from  "read_parsers.hh":
    cdef cppclass CpSequence "khmer::read_parsers::Read":
        string name
        string annotations
        string sequence
        string quality

        void reset()

    ctypedef pair[CpSequence,CpSequence] CpSequencePair "khmer::read_parsers::ReadPair"

    cdef cppclass CpFastxParser "khmer::read_parsers::FastxParser":
        CpFastxParser(const char *)
        bool is_complete()
        void imprint_next_read(CpSequence&) except +
        void imprint_next_read_pair(CpSequencePair&, uint8_t) except +



########################################################################
#
# Hashtable: Bindings for the existing CPython Hashtable wrapper.
#
########################################################################

# All we really need are the PyObject struct definitions
# for our extension objects.
cdef extern from "_khmer.hh":
    ctypedef struct CPyHashtable_Object "khmer::khmer_KHashtable_Object":
        CpHashtable* hashtable


cdef extern from "hashtable.hh" namespace "khmer":
    cdef cppclass CpHashtable "khmer::Hashtable":
        CpHashtable(WordLength)
        const BoundedCounterType get_count(const char *) const
        const BoundedCounterType get_count(HashIntoType) const
        const WordLength ksize() const


########################################################################
#
# Traversal: wrapper for traversal.hh.
#
########################################################################


cdef extern from "traversal.hh":
    cdef cppclass CpTraverser "khmer::Traverser":
        CpTraverser(CpHashtable *)

        void push_filter(KmerFilter)
        KmerFilter pop_filter()
    
        uint32_t traverse(const CpKmer&, KmerQueue&) const
        uint32_t traverse_left(const CpKmer&, KmerQueue&) const     
        uint32_t traverse_right(const CpKmer&, KmerQueue&) const

        uint32_t degree(const CpKmer&) const
        uint32_t degree_left(const CpKmer&) const
        uint32_t degree_right(const CpKmer&) const


########################################################################
#
# Partitioning: Component and StreamingPartitioner
#
########################################################################


cdef extern from "partitioning.hh" namespace "khmer":

    cdef cppclass CpComponent "khmer::Component":
        CpComponent()
        CpComponent(uint64_t)

        const uint64_t component_id
        set[HashIntoType] tags

        void merge(set[shared_ptr[CpComponent]])
        void add_tag(HashIntoType)
        void add_tag(set[HashIntoType])
        uint64_t get_n_tags() const
        uint64_t get_n_created() const
        uint64_t get_n_destroyed() const

    ctypedef shared_ptr[CpComponent] ComponentPtr
    ctypedef set[ComponentPtr] ComponentPtrSet

    cdef cppclass CpGuardedKmerCompMap "khmer::GuardedKmerCompMap":
        map[HashIntoType, ComponentPtr] data

        ComponentPtr get(HashIntoType)
        void set(HashIntoType, ComponentPtr)
        bool contains(HashIntoType)

    cdef cppclass CpStreamingPartitioner "khmer::StreamingPartitioner":
        CpStreamingPartitioner(CpHashtable * ) except +MemoryError

        void consume(string&) except +MemoryError
        void consume_pair(string&, string&) except +MemoryError
        uint64_t consume_fasta(string&) except +MemoryError

        void add_component(ComponentPtr comp)
        void map_tags_to_component(set[HashIntoType]&, ComponentPtr&)
        void find_connected_tags(queue[CpKmer]&, 
                                 set[HashIntoType]&,
                                 set[HashIntoType]&) except +MemoryError
        uint64_t get_n_components() const
        uint64_t get_n_tags() const

        ComponentPtr get_tag_component(string&) const
        ComponentPtr get_nearest_component(string&) const

        weak_ptr[ComponentPtrSet] get_component_set()
        weak_ptr[CpGuardedKmerCompMap] get_tag_component_map()


