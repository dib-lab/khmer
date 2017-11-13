from libcpp cimport bool
from libcpp.memory cimport shared_ptr
from libcpp.queue cimport queue
from libcpp.set cimport set
from libcpp.string cimport string

from khmer._oxli.oxli_types cimport *

cdef extern from "oxli/kmer_hash.hh" namespace "oxli":
    cdef cppclass CpKmer "oxli::Kmer":
        HashIntoType kmer_f
        HashIntoType kmer_r
        HashIntoType kmer_u

        CpKmer(HashIntoType, HashIntoType, HashIntoType)
        CpKmer(string, WordLength)
        CpKmer(const CpKmer&)
        CpKmer()

        bool is_forward() const
        void set_from_unique_hash(HashIntoType, WordLength)

    cdef cppclass CpKmerFactory "oxli::KmerFactory":
        KmerFactory(WordLength)

        CpKmer build_kmer(HashIntoType) const
        CpKmer build_kmer(HashIntoType, HashIntoType) const
        CpKmer build_kmer(string &) const
        CpKmer build_kmer(const char *) const

    cdef cppclass CpKmerIterator "oxli::KmerIterator" (CpKmerFactory):
        CpKmerIterator(const char *, unsigned char)
        CpKmer first(HashIntoType &, HashIntoType &)
        CpKmer next(HashIntoType &, HashIntoType &)
        CpKmer first()
        CpKmer next()
        bool done()
        unsigned int get_start_pos() const
        unsigned int get_end_pos() const

    HashIntoType _hash(const string, const WordLength)
    HashIntoType _hash(const string, const WordLength, 
                       HashIntoType &, HashIntoType &)
    HashIntoType _hash(const char *, const WordLength)
    HashIntoType _hash(const char *, const WordLength,
                       HashIntoType &, HashIntoType &)
    HashIntoType _hash_forward(const char *, WordLength)
    string _revhash(HashIntoType, WordLength)
    string _revcomp(const string&)
    HashIntoType _hash_murmur(const string&, const WordLength)
    HashIntoType _hash_murmur(const string&,
                              HashIntoType&, HashIntoType&)
    HashIntoType _hash_murmur_forward(const string&)


cdef extern from "oxli/oxli.hh" namespace "oxli":
    ctypedef queue[CpKmer] KmerQueue
    ctypedef set[CpKmer] KmerSet
    ctypedef bool (*KmerFilter) (CpKmer kmer)


cdef class Kmer:
    cdef shared_ptr[CpKmer] _this
    cdef readonly str kmer

    @staticmethod
    cdef Kmer wrap(CpKmer * cpkmer, WordLength K)
