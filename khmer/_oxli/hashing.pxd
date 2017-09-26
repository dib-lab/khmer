from libcpp cimport bool
from libcpp.memory cimport shared_ptr
from libcpp.queue cimport queue
from libcpp.set cimport set
from libcpp.string cimport string
from libc.stdint cimport uint64_t

from khmer._oxli.oxli_types cimport *
from khmer._oxli.utils cimport oxli_raise_py_error

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


cdef extern from "oxli/hashing.hh" namespace "oxli":
    cdef cppclass CpHashedKmer "oxli::HashedKmer":
        HashIntoType kmer_fw
        HashIntoType kmer_rc

        CpHashedKmer(HashIntoType, HashIntoType)
        CpHashedKmer()

        bool forward() const
        bool operator< (const CpHashedKmer&) const
        HashIntoType key() const

    cdef cppclass HashedKmerIterator:
        CpHashedKmer next() except +oxli_raise_py_error
        bool done() const

    cdef cppclass KmerHasher

    cdef cppclass CpSequenceHasher "oxli::SequenceHasher" (HashedKmerIterator):
        CpSequenceHasher(string, const KmerHasher *)
        uint64_t get_start_pos() const
        uint64_t get_end_pos() const

    cdef cppclass KmerHasher:
        KmerHasher()
        
        CpHashedKmer hash_dna(const char *) const
        CpHashedKmer hash_dna(const string&) const
        CpHashedKmer hash_dna_fw(const string&) const
        CpHashedKmer hash_dna_fw(const char *) const
        CpSequenceHasher hash_sequence(const string&) const

    cdef cppclass ReversibleKmerHasher (KmerHasher):
        ReversibleKmerHasher()

        string unhash_dna(CpHashedKmer& ) const

    cdef cppclass ShiftingKmerHasher (KmerHasher):
        ShiftingKmerHasher()

        CpHashedKmer _get_left(const CpHashedKmer&, const char) const
        CpHashedKmer _right_right(const CpHashedKmer&, const char) const

    cdef cppclass CpTwoBitKmerHasher "oxli::TwoBitKmerHasher" \
        (ReversibleKmerHasher, ShiftingKmerHasher):

        CpTwoBitKmerHasher(WordLength)


cdef class SequenceHasher:
    cdef shared_ptr[CpSequenceHasher] _sh_this
    cdef shared_ptr[CpTwoBitKmerHasher] _hasher
    cdef string _sequence
    cdef readonly WordLength K


cdef extern from "oxli/oxli.hh" namespace "oxli":
    ctypedef queue[CpKmer] KmerQueue
    ctypedef set[CpKmer] KmerSet
    ctypedef bool (*KmerFilter) (CpKmer kmer)


cdef class Kmer:
    cdef shared_ptr[CpKmer] _this
    cdef readonly str kmer

    @staticmethod
    cdef Kmer wrap(CpKmer * cpkmer, WordLength K)
