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

####

