from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.set cimport set
from libc.stdint cimport uint32_t, uint8_t


cdef extern from "khmer.hh" namespace "khmer":
  ctypedef unsigned long long int ExactCounterType
  ctypedef unsigned long long int HashIntoType
  ctypedef unsigned char WordLength
  ctypedef unsigned short int BoundedCounterType
  ctypedef void (*CallbackFn)(const char *, void *,
                              unsigned long long,
                              unsigned long long)
  #ctypedef set[HashIntoType] SeenSet
  ctypedef set[unsigned long long int] SeenSet
  ctypedef unsigned int PartitionID


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
        SubsetPartition * partition
        unsigned int consume_string(const string &)


cdef extern from "counting.hh" namespace "khmer":
    cdef cppclass CountingHash:
        CountingHash(WordLength, HashIntoType, uint32_t)
        CountingHash(WordLength, vector[unsigned long long int]&, uint32_t)
        void get_kadian_count(const string &, BoundedCounterType &, unsigned int)
        HashIntoType * fasta_count_kmers_by_position(string &,
                        const unsigned int,
                        BoundedCounterType,
                        CallbackFn,
                        void *)
        BoundedCounterType get_max_count(const string &s)
        BoundedCounterType get_min_count(const string &s)
        void output_fasta_kmer_pos_freq(const string &, const string &)

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
        unsigned int consume_high_abund_kmers(const string &, BoundedCounterType)
        const HashIntoType n_occupied(HashIntoType, HashIntoType)

        HashIntoType * abundance_distribution(string, CppHashbits *)
        unsigned int trim_on_abundance(string, BoundedCounterType) const
        unsigned int trim_below_abundance(string, BoundedCounterType) const
        void set_use_bigcount(bool)
        void count(const char *)
        void count(HashIntoType)
        vector[HashIntoType] get_tablesizes() const


cdef extern from "hashbits.hh" namespace "khmer":
    cdef cppclass CppHashbits "khmer::Hashbits":
        CppHashbits(WordLength, vector[unsigned long long int]&)
        const HashIntoType n_kmers(HashIntoType, HashIntoType) const

        # FIXME: this is from hashtable, how to avoid redeclaring
        # all inherited methods?
        SubsetPartition * partition
        SeenSet all_tags
        const WordLength ksize() const
        vector[HashIntoType] get_tablesizes() const
        void get_median_count(const string &, BoundedCounterType &, float &, float &)
        void save(string)
        void consume_fasta(const string &,
                           unsigned int &,
                           unsigned long long &,
                           CallbackFn,
                           void *)
        void consume_fasta_and_tag(string &,
                                   unsigned int &,
                                   unsigned long long &,
                                   CallbackFn,
                                   void *)
        void consume_fasta_and_tag_with_stoptags(const string &,
                                                 unsigned int &,
                                                 unsigned long long &,
                                                 CallbackFn,
                                                 void *)
        void calc_connected_graph_size(const char *,
                                       unsigned long long&,
                                       SeenSet&,
                                       const unsigned long long,
                                       bool) const
        unsigned int consume_string(const string &)
        void add_kmer_to_tags(HashIntoType)
        unsigned int _get_tag_density() const
        void _set_tag_density(unsigned int)
        const HashIntoType n_occupied(HashIntoType, HashIntoType)
        const BoundedCounterType get_count(const char *) const
        const BoundedCounterType get_count(HashIntoType) const
        void count(const char *)
        void count(HashIntoType)
        void filter_if_present(const string,
                               const string,
                               CallbackFn,
                               void *)
        void consume_partitioned_fasta(const string &,
                                       unsigned int &,
                                       unsigned long long &,
                                       CallbackFn,
                                       void *)
        unsigned int count_kmers_within_radius(HashIntoType,
                                               HashIntoType,
                                               unsigned int,
                                               unsigned int,
                                               const SeenSet *) const
        unsigned int count_kmers_on_radius(HashIntoType,
                                           HashIntoType,
                                           unsigned int,
                                           unsigned int) const
        unsigned int kmer_degree(HashIntoType, HashIntoType) const
        unsigned int kmer_degree(const char *) const
        unsigned int find_radius_for_volume(HashIntoType,
                                            HashIntoType,
                                            unsigned int,
                                            unsigned int) const
        void add_tag(HashIntoType)
        void add_stop_tag(HashIntoType)
        void save_tagset(string)
        void load_tagset(string, bool)
        void identify_stop_tags_by_position(string, vector[unsigned int] &) const
        void extract_unique_paths(string, unsigned int, float, vector[string]&)


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


cdef extern from "subset.hh" namespace "khmer":
    cdef cppclass SubsetPartition:
        SubsetPartition(CppHashbits *)
        void do_partition(HashIntoType, HashIntoType,
                          bool, bool,
                          CallbackFn, void *)
        void count_partitions(unsigned int&, unsigned int&)
        unsigned int output_partitioned_file(const string,
                                             const string,
                                             bool,
                                             CallbackFn,
                                             void *)
        void merge(SubsetPartition *)
        void find_all_tags(HashIntoType, HashIntoType,
                           SeenSet&, const SeenSet&, bool, bool)
        PartitionID assign_partition_id(HashIntoType, SeenSet&)
        PartitionID get_partition_id(string)
        PartitionID get_partition_id(HashIntoType)
        PartitionID join_partitions(PartitionID, PartitionID)
        unsigned int find_unpart(const string, bool, bool, CallbackFn, void *)

    cdef cppclass pre_partition_info:
        HashIntoType kmer
        SeenSet tagged_kmers

        pre_partition_info(HashIntoType)
        void count_partitions(unsigned int&, unsigned int&)
