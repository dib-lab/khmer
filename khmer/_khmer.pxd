from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.set cimport set
from libc.stdint cimport uint32_t, uint8_t, uint64_t


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

  ctypedef unsigned long long int Label
  ctypedef map[Label, Label*] LabelPtrMap
  ctypedef set[Label*] LabelPtrSet


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
    cdef cppclass CppHashtable "khmer::Hashtable":
        CppHashtable(WordLength, uint32_t const, uint8_t const)
        SubsetPartition * partition
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

        void count(const char *)
        void count(HashIntoType)

        SubsetPartition * partition
        SeenSet all_tags
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
        void add_kmer_to_tags(HashIntoType)
        unsigned int _get_tag_density() const
        void _set_tag_density(unsigned int)
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
        HashIntoType * abundance_distribution(string, CppHashbits *)
        unsigned int trim_on_abundance(string, BoundedCounterType) const
        unsigned int trim_below_abundance(string, BoundedCounterType) const
        void set_use_bigcount(bool)
        vector[HashIntoType] get_tablesizes() const
        const HashIntoType n_occupied(HashIntoType, HashIntoType)


cdef extern from "hashbits.hh" namespace "khmer":
    cdef cppclass CppHashbits "khmer::Hashbits":
        CppHashbits(WordLength, vector[unsigned long long int]&)
        const HashIntoType n_kmers(HashIntoType, HashIntoType) const
        vector[HashIntoType] get_tablesizes() const
        const HashIntoType n_occupied(HashIntoType, HashIntoType)


cdef extern from "labelhash.hh" namespace "khmer":
    cdef cppclass CppLabelHash "khmer::LabelHash":
        CppLabelHash(WordLength, vector[unsigned long long int]&)
        LabelPtrMap label_ptrs
        unsigned int n_labels() const
        void consume_fasta_and_tag_with_labels(string &filename,
                                               unsigned int &,
                                               unsigned long long &,
                                               CallbackFn,
                                               void *)
        void consume_partitioned_fasta_and_tag_with_labels(string &filename,
                                                           unsigned int &,
                                                           unsigned long long &,
                                                           CallbackFn,
                                                           void *)
        LabelPtrSet get_tag_labels(const HashIntoType&)
        unsigned int sweep_label_neighborhood(const string &,
                                              LabelPtrSet&,
                                              unsigned int,
                                              bool,
                                              bool)
        void consume_sequence_and_tag_with_labels(const string&,
                                                  unsigned long long&,
                                                  Label&,
                                                  SeenSet *)
        Label * check_and_allocate_label(Label)

        # FIXME: this is from hashtable, how to avoid redeclaring
        # all inherited methods?
        SubsetPartition * partition
        SeenSet all_tags
        const WordLength ksize() const
        const BoundedCounterType get_count(const char *) const
        const BoundedCounterType get_count(HashIntoType) const
        const HashIntoType n_occupied(HashIntoType, HashIntoType)
        unsigned int consume_string(const string &)
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
        unsigned int _get_tag_density() const
        void _set_tag_density(unsigned int)
        void count(const char *)
        void count(HashIntoType)

        # FIXME: this is from hashbits, how to avoid redeclaring
        # all inherited methods?
        const HashIntoType n_kmers(HashIntoType, HashIntoType) const


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
        SubsetPartition(CppHashtable *)
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
        unsigned int sweep_for_tags(const string&,
                        SeenSet& tagged_kmers,
                        const SeenSet& all_tags,
                        unsigned int range,
                        bool break_on_stop_tags,
                        bool stop_big_traversals)

    cdef cppclass pre_partition_info:
        HashIntoType kmer
        SeenSet tagged_kmers

        pre_partition_info(HashIntoType)
        void count_partitions(unsigned int&, unsigned int&)


cdef extern from "khmer_config.hh" namespace "khmer":
    cdef cppclass CppConfig "khmer::Config":
        CppConfig()
        bool has_extra_sanity_checks() const

        uint32_t get_number_of_threads() const
        void set_number_of_threads(uint32_t)

        uint64_t get_reads_input_buffer_size() const
        void set_reads_input_buffer_size(uint64_t)

        uint8_t get_input_buffer_trace_level() const
        void set_input_buffer_trace_level(uint8_t const)
        uint8_t get_reads_parser_trace_level() const
        void set_reads_parser_trace_level(uint8_t const)
    cdef CppConfig & get_active_config()
