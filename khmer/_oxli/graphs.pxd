from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.set cimport set
from libcpp.memory cimport unique_ptr
from libc.stdint cimport uint8_t, uint32_t, uint64_t, uintptr_t

from oxli_types cimport *
from hashing cimport CpKmer, KmerSet
from parsing cimport CpReadParser, CpSequence

# All we really need are the PyObject struct definitions
# for our extension objects.
cdef extern from "khmer/_cpy_khmer.hh":

    ctypedef struct CPyHashtable_Object "khmer::khmer_KHashtable_Object":
        CpHashtable * hashtable

    ctypedef struct CPyHashgraph_Object "khmer::khmer_KHashgraph_Object":
        CPyHashtable_Object khashtable
        CpHashgraph * hashgraph

    ctypedef struct CPyNodegraph_Object "khmer::khmer_KNodegraph_Object":
        CPyHashgraph_Object khashgraph
        CpNodegraph * nodegraph

    ctypedef struct CPyCountgraph_Object "khmer::khmer_KCountgraph_Object":
        CPyHashgraph_Object khashgraph
        CpCountgraph * countgraph

    ctypedef struct CPyGraphLabels_Object "khmer::khmer_KGraphLabels_Object":
        CpLabelHash * labelhash


cdef extern from "oxli/hashtable.hh" namespace "oxli":
    cdef cppclass CpHashtable "oxli::Hashtable":
        const WordLength ksize() const
        HashIntoType hash_dna(const char *) const
        HashIntoType hash_dna_top_strand(const char *) const
        HashIntoType hash_dna_bottom_strand(const char *) const
        string unhash_dna(HashIntoType) const
        void count(const char *)
        void count(HashIntoType)
        bool add(const char *)
        bool add(HashIntoType)
        const BoundedCounterType get_count(const char *) const
        const BoundedCounterType get_count(HashIntoType) const
        void save(string)
        void load(string)
        uint32_t consume_string(const string &)
        bool check_and_normalize_read(string &) const
        uint32_t check_and_process_read(string &, bool &)
        void consume_seqfile[SeqIO](const string &, uint32_t &, uint64_t &)
        void consume_seqfile[SeqIO](unique_ptr[CpReadParser[SeqIO]]&,
                                    uint32_t &, uint64_t &)
        void set_use_bigcount(bool)
        bool get_use_bigcount()
        bool median_at_least(const string &, uint32_t cutoff)
        void get_median_count(const string &, BoundedCounterType &,
                              float &, float &)
        const uint64_t n_unique_kmers() const
        const uint64_t n_occupied() const
        vector[uint64_t] get_tablesizes() const
        const uintptr_t n_tables() const
        void get_kmers(const string &, vector[string] &)
        void get_kmer_hashes(const string &, vector[HashIntoType] &) const
        void get_kmer_hashes_as_hashset(const string &,
                                        set[HashIntoType]) const
        void get_kmer_counts(const string &,
                             vector[BoundedCounterType] &) const
        uint8_t ** get_raw_tables()
        BoundedCounterType get_min_count(const string &)
        BoundedCounterType get_max_count(const string &)
        uint64_t * abundance_distribution[SeqIO](unique_ptr[CpReadParser[SeqIO]]&,
                                          CpHashtable *)
        uint64_t * abundance_distribution[SeqIO](string, CpHashtable *)
        uint64_t trim_on_abundance(string, BoundedCounterType) const
        uint64_t trim_below_abundance(string, BoundedCounterType) const
        vector[uint32_t] find_spectral_error_positions(string,
                                                       BoundedCounterType)

    cdef cppclass CpCounttable "oxli::Counttable" (CpHashtable):
        CpCounttable(WordLength, vector[uint64_t])

    cdef cppclass CpNodetable "oxli::Nodetable" (CpHashtable):
        CpNodetable(WordLength, vector[uint64_t])

    cdef cppclass CpQFCounttable "oxli::QFCounttable":
        CpQFCounttable(WordLength, int)
        void count(const char *)
        void count(HashIntoType)
        bool add(const char *)
        bool add(HashIntoType)
        const BoundedCounterType get_count(const char *) const
        const BoundedCounterType get_count(HashIntoType) const
        HashIntoType hash_dna(const char *) const
        string unhash_dna(HashIntoType) const
        const WordLength ksize() const
        vector[uint64_t] get_tablesizes() const
        void get_kmers(const string &, vector[string] &)
        uint32_t consume_string(const string &)
        void get_kmer_counts(const string &,
                             vector[BoundedCounterType] &) const
        void get_kmer_hashes(const string &, vector[HashIntoType] &) const
        BoundedCounterType get_min_count(const string &)
        BoundedCounterType get_max_count(const string &)
        void get_median_count(const string &, BoundedCounterType &,
                              float &, float &)
        uint64_t * abundance_distribution[SeqIO](string, CpHashtable *)
        uint64_t * abundance_distribution[SeqIO](unique_ptr[CpReadParser[SeqIO]]&,
                                                 CpHashtable *)
        uint64_t trim_on_abundance(string, BoundedCounterType) const
        uint64_t trim_below_abundance(string, BoundedCounterType) const
        vector[uint32_t] find_spectral_error_positions(string,
                                                       BoundedCounterType)
        void consume_seqfile[SeqIO](const string &, uint32_t &, uint64_t &)
        void consume_seqfile[SeqIO](unique_ptr[CpReadParser[SeqIO]]&,
                                    uint32_t &, uint64_t &)
        void save(string)
        void load(string)
        const uint64_t n_unique_kmers() const
        const uint64_t n_occupied() const
        const uintptr_t n_tables() const


cdef extern from "oxli/hashgraph.hh" namespace "oxli":
    cdef cppclass CpHashgraph "oxli::Hashgraph" (CpHashtable):
        void _set_tag_density(unsigned int)
        unsigned int _get_tag_density() const
        void add_tag(HashIntoType)
        void add_stop_tag(HashIntoType)
        uintptr_t n_tags() const
        void divide_tags_into_subsets(unsigned int, set[HashIntoType] &)
        void add_kmer_to_tags(HashIntoType)
        void clear_tags()
        void consume_seqfile_and_tag[SeqIO](unique_ptr[CpReadParser[SeqIO]]&,
                                   unsigned int &,
                                   unsigned long long)
        void consume_seqfile_and_tag[SeqIO](const string &,
                                   unsigned int &,
                                   unsigned long long &)
        void consume_sequence_and_tag(const string &,
                                      unsigned long long &,
                                      set[HashIntoType] &)
        void consume_partitioned_fasta[SeqIO](const string &,
                                       unsigned int &,
                                       unsigned long long &)
        uintptr_t trim_on_stoptags(string) const
        unsigned int traverse_from_kmer(CpKmer,
                                        uint32_t,
                                        KmerSet&,
                                        uint32_t) const
        void print_tagset(string)
        void save_tagset(string)
        void load_tagset(string)
        void print_stop_tags(string)
        void save_stop_tags(string)
        void load_stop_tags(string)
        void load_stop_tags(string, bool)
        void extract_unique_paths(string, uint32_t, float, vector[string])
        void calc_connected_graph_size(CpKmer, uint64_t&, KmerSet&,
                                       const uint64_t, bool) const
        uint32_t kmer_degree(HashIntoType, HashIntoType)
        uint32_t kmer_degree(const char *)
        void find_high_degree_nodes(const char *, set[HashIntoType] &) const
        unsigned int traverse_linear_path(const CpKmer,
                                          set[HashIntoType] &,
                                          set[HashIntoType] &,
                                          CpHashtable &,
                                          set[HashIntoType] &) const
        void _validate_pmap()

    cdef cppclass CpCountgraph "oxli::Countgraph" (CpHashgraph):
        CpCountgraph(WordLength, vector[uint64_t])

    cdef cppclass CpNodegraph "oxli::Nodegraph" (CpHashgraph):
        CpNodegraph(WordLength, vector[uint64_t])


cdef extern from "oxli/labelhash.hh" namespace "oxli":
    cdef cppclass CpLabelHash "oxli::LabelHash":
        CpLabelHash(CpHashgraph *)
        uintptr_t n_labels() const
        void consume_seqfile_and_tag_with_labels[SeqIO](const string &,
                                               uint32_t &,
                                               uint64_t &,
                                               CallbackFn,
                                               void *)
        void consume_seqfile_and_tag_with_labels[SeqIO](const string &,
                                               uint32_t &,
                                               uint64_t &)
        void consume_seqfile_and_tag_with_labels[SeqIO](
                               unique_ptr[CpReadParser[SeqIO]]&,
                               uint32_t &,
                               uint64_t &,
                               CallbackFn,
                               void *)
        void consume_seqfile_and_tag_with_labels[SeqIO](
                               unique_ptr[CpReadParser[SeqIO]]&,
                               uint32_t &,
                               uint64_t &)
        void get_tag_labels(const HashIntoType, LabelSet &)
        void get_tags_from_label(const Label,
                                  TagSet&)
        void link_tag_and_label(const HashIntoType, const Label)
        uint32_t sweep_label_neighborhood(const string &,
                                           LabelSet &,
                                           uint32_t,
                                           bool,
                                           bool)
        void save_labels_and_tags(string)
        void load_labels_and_tags(string)
        void label_across_high_degree_nodes(const char *,
                                             set[HashIntoType] &,
                                             const Label)

cdef CpHashgraph * get_hashgraph_ptr(object graph)
cdef CpLabelHash * get_labelhash_ptr(object graph)
