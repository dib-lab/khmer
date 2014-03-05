from libcpp.string cimport string


cdef extern from "khmer.hh" namespace "khmer":
  ctypedef unsigned long long int ExactCounterType
  ctypedef unsigned long long int HashIntoType
  ctypedef unsigned char WordLength


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
