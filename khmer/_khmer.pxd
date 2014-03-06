from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map
from libc.stdint cimport uint32_t, uint8_t


cdef extern from "khmer.hh" namespace "khmer":
  ctypedef unsigned long long int ExactCounterType
  ctypedef unsigned long long int HashIntoType
  ctypedef unsigned char WordLength
  ctypedef unsigned short int BoundedCounterType
  ctypedef void (*CallbackFn)(const char *, void *,
                              unsigned long long,
                              unsigned long long)


cdef extern from "ktable.hh" namespace "khmer":
    cdef cppclass CppKTable "khmer::KTable":
        CppKTable(long)
        ExactCounterType get_count(const char *)
        ExactCounterType get_count(HashIntoType)
        void count(const char *)
        void set_count(const char *, ExactCounterType c)
        void set_count(HashIntoType, ExactCounterType c)

        HashIntoType n_entries()
        const WordLength ksize() const
        const HashIntoType max_hash() const

        void consume_string(const string &)
        void clear()
        void update(const CppKTable &)
        CppKTable * intersect(const CppKTable &) const

    cdef HashIntoType _hash(const char*, const WordLength)
    cdef HashIntoType _hash(const char *, const WordLength,
                            HashIntoType&, HashIntoType&)
    cdef HashIntoType _hash_forward(const char *, WordLength)
    cdef string _revhash(HashIntoType, WordLength)


cdef extern from "hashtable.hh" namespace "khmer":
    cdef cppclass Hashtable:
        Hashtable(WordLength, uint32_t const, uint8_t const)
        unsigned int consume_string(const string &)


cdef extern from "counting.hh" namespace "khmer":
    cdef cppclass CountingHash:
        CountingHash(WordLength, HashIntoType, uint32_t)
        CountingHash(WordLength, vector[unsigned long long int]&, uint32_t)
        void get_kadian_count(const string &, BoundedCounterType &, unsigned int)

        # FIXME: this is from hashtable, how to avoid redeclaring
        # all inherited methods?
        const WordLength ksize() const
        unsigned int consume_string(const string &)
        const BoundedCounterType get_count(const char *) const
        const BoundedCounterType get_count(HashIntoType) const
        void get_median_count(const string &, BoundedCounterType &, float &, float &)
        void consume_fasta(const string &,
                           unsigned int &,
                           unsigned long long &,
                           CallbackFn,
                           void *)
        void save(string)
        void load(string)


cdef extern from "aligner.hh" namespace "khmer":
    cdef cppclass CandidateAlignment:
        CandidateAlignment(map[int,int], string)
        CandidateAlignment()
        string getReadAlignment(string)
        string alignment

    cdef cppclass Aligner:
        Aligner(CountingHash*, double, double, unsigned int)
        CandidateAlignment align(CountingHash*, const string&,
                                 const string&, int)
        CandidateAlignment alignRead(const string&)
