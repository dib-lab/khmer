from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.set cimport set
from libcpp.memory cimport unique_ptr, shared_ptr, weak_ptr
from libc.stdint cimport uint8_t, uint32_t, uint64_t, uintptr_t

from khmer._oxli.oxli_types cimport *
from khmer._oxli.hashing cimport Kmer, CpKmer, KmerSet, CpKmerFactory, CpKmerIterator
from khmer._oxli.parsing cimport CpReadParser, CpSequence, FastxParserPtr
from khmer._oxli.legacy_partitioning cimport (CpSubsetPartition, cp_pre_partition_info,
                                   SubsetPartition)
from khmer._oxli.utils cimport oxli_raise_py_error


cdef extern from "Python.h":
    ctypedef struct PyObject
    object PyMemoryView_FromBuffer(Py_buffer *view)


cdef extern from "oxli/storage.hh":
    cdef cppclass CpStorage "oxli::Storage":
        CpStorage()

        vector[uint64_t] get_tablesizes()
        const size_t n_tables()
        void save(string, WordLength)
        void load(string, WordLength&)
        const uint64_t n_occupied()
        const uint64_t n_unique_kmers()
        BoundedCounterType test_and_set_bits(HashIntoType)
        bool add(HashIntoType khash)
        const BoundedCounterType get_count(HashIntoType)
        uint8_t ** get_raw_tables()

        void set_use_bigcount(bool)
        bool get_use_bigcount()


cdef extern from "oxli/hashtable.hh" namespace "oxli" nogil:
    cdef cppclass CpHashtable "oxli::Hashtable" (CpKmerFactory):
        const WordLength ksize() const
        HashIntoType hash_dna(const char *) except +oxli_raise_py_error
        HashIntoType hash_dna_top_strand(const char *) except +oxli_raise_py_error
        HashIntoType hash_dna_bottom_strand(const char *) except +oxli_raise_py_error
        string unhash_dna(HashIntoType) except +oxli_raise_py_error
        void count(const char *)
        void count(HashIntoType)
        bool add(const char *)
        bool add(HashIntoType)
        const BoundedCounterType get_count(const char *) except +oxli_raise_py_error
        const BoundedCounterType get_count(HashIntoType) except +oxli_raise_py_error
        void save(string)
        void load(string) except +oxli_raise_py_error
        uint32_t consume_string(const string &)
        bool check_and_normalize_read(string &) const
        uint32_t check_and_process_read(string &, bool &)

        void consume_seqfile[SeqIO](const string &, uint32_t &, uint64_t &) except +oxli_raise_py_error
        void consume_seqfile[SeqIO](shared_ptr[CpReadParser[SeqIO]]&,
                                    uint32_t &, uint64_t &) except +oxli_raise_py_error

        void consume_seqfile_with_mask[SeqIO](const string &, CpHashtable *,
                                              uint32_t, uint32_t &, uint64_t &, bool) except +oxli_raise_py_error
        void consume_seqfile_with_mask[SeqIO](shared_ptr[CpReadParser[SeqIO]]&, CpHashtable *,
                                              uint32_t, uint32_t &, uint64_t &, bool) except +oxli_raise_py_error

        void consume_seqfile_banding[SeqIO](const string &, uint32_t, uint32_t, uint32_t &, uint64_t &) except +oxli_raise_py_error
        void consume_seqfile_banding[SeqIO](shared_ptr[CpReadParser[SeqIO]]&,
                                    uint32_t, uint32_t, uint32_t &, uint64_t &) except +oxli_raise_py_error

        void consume_seqfile_banding_with_mask[SeqIO](const string &, uint32_t, uint32_t,
                                                      CpHashtable *, uint32_t, uint32_t &,
                                                      uint64_t &, bool) except +oxli_raise_py_error
        void consume_seqfile_banding_with_mask[SeqIO](shared_ptr[CpReadParser[SeqIO]]&,
                                                      uint32_t, uint32_t,
                                                      CpHashtable *, uint32_t,
                                                      uint32_t &, uint64_t &, bool) except +oxli_raise_py_error

        void set_use_bigcount(bool) except +ValueError
        bool get_use_bigcount()
        bool median_at_least(const string &, uint32_t cutoff)
        void get_median_count(const string &, BoundedCounterType &,
                              float &, float &) except +oxli_raise_py_error
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
        uint64_t * abundance_distribution[SeqIO](string, CpHashtable *) except +oxli_raise_py_error
        uint64_t * abundance_distribution[SeqIO](shared_ptr[CpReadParser[SeqIO]]&,
                                          CpHashtable *) except +oxli_raise_py_error
        uint64_t trim_on_abundance(string, BoundedCounterType) const
        uint64_t trim_below_abundance(string, BoundedCounterType) const
        vector[uint32_t] find_spectral_error_positions(string,
                                                       BoundedCounterType)

    cdef cppclass CpMurmurHashtable "oxli::MurmurHashtable" (CpHashtable):
        CpMurmurHashtable(WordLength, CpStorage *)

    cdef cppclass CpCyclicHashtable "oxli::CyclicHashtable" (CpHashtable):
        CpCyclicHashtable(WordLength, CpStorage *)

    cdef cppclass CpCounttable "oxli::Counttable" (CpMurmurHashtable):
        CpCounttable(WordLength, vector[uint64_t])

    cdef cppclass CpCyclicCounttable "oxli::CyclicCounttable" (CpCyclicHashtable):
        CpCyclicCounttable(WordLength, vector[uint64_t])

    cdef cppclass CpSmallCounttable "oxli::SmallCounttable" (CpMurmurHashtable):
        CpSmallCounttable(WordLength, vector[uint64_t])

    cdef cppclass CpNodetable "oxli::Nodetable" (CpMurmurHashtable):
        CpNodetable(WordLength, vector[uint64_t])

    cdef cppclass CpQFCounttable "oxli::QFCounttable" (CpHashtable):
        CpQFCounttable(WordLength, uint64_t) except +oxli_raise_py_error


cdef extern from "oxli/hashgraph.hh" namespace "oxli" nogil:
    cdef cppclass CpHashgraph "oxli::Hashgraph" (CpHashtable):
        set[HashIntoType] all_tags
        set[HashIntoType] stop_tags
        set[HashIntoType] repart_small_tags
        shared_ptr[CpSubsetPartition] partition

        void _set_tag_density(unsigned int)
        unsigned int _get_tag_density() const
        void add_tag(HashIntoType)
        void add_stop_tag(HashIntoType)
        bool has_tag(HashIntoType)
        bool hash_stop_tag(HashIntoType)
        uintptr_t n_tags() const
        void divide_tags_into_subsets(unsigned int, set[HashIntoType] &)
        void add_kmer_to_tags(HashIntoType) nogil
        void clear_tags()

        void consume_seqfile_and_tag[SeqIO](const string &,
                                   unsigned int,
                                   unsigned long long)

        # Ugly workaround. For some reason, Cython doesn't like *just this*
        # templated overload -- it chooses whichever was defined last, breaking
        # resolution for either strings of FastxParserPtr. So, we rename it on
        # the Cython side and give it a real name substitution for code gen.
        void consume_seqfile_and_tag_readparser "consume_seqfile_and_tag" [SeqIO](shared_ptr[CpReadParser[SeqIO]],
                                   unsigned int,
                                   unsigned long long)

        void consume_sequence_and_tag(const string &,
                                      unsigned long long &)

        void consume_sequence_and_tag(const string &,
                                      unsigned long long &,
                                      set[HashIntoType] &)

        void consume_partitioned_fasta[SeqIO](const string &,
                                       unsigned int &,
                                       unsigned long long &) except +oxli_raise_py_error

        uintptr_t trim_on_stoptags(string)

        unsigned int traverse_from_kmer(CpKmer,
                                        uint32_t,
                                        KmerSet&,
                                        uint32_t) nogil
        void get_tags_for_sequence(string&, set[HashIntoType]&)
        void print_tagset(string)
        void save_tagset(string)
        void load_tagset(string) except +oxli_raise_py_error
        void load_tagset(string, bool) except +oxli_raise_py_error
        void print_stop_tags(string)
        void save_stop_tags(string)
        void load_stop_tags(string) except +oxli_raise_py_error
        void load_stop_tags(string, bool) except +oxli_raise_py_error
        void extract_unique_paths(string, uint32_t, float, vector[string])
        void calc_connected_graph_size(CpKmer, uint64_t&, KmerSet&,
                                       const uint64_t, bool)
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

    cdef cppclass CpSmallCountgraph "oxli::SmallCountgraph" (CpHashgraph):
        CpSmallCountgraph(WordLength, vector[uint64_t])

    cdef cppclass CpNodegraph "oxli::Nodegraph" (CpHashgraph):
        CpNodegraph(WordLength, vector[uint64_t])

        void update_from(const CpNodegraph &) except +oxli_raise_py_error


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
                               shared_ptr[CpReadParser[SeqIO]]&,
                               uint32_t &,
                               uint64_t &,
                               CallbackFn,
                               void *)
        void consume_seqfile_and_tag_with_labels[SeqIO](
                               shared_ptr[CpReadParser[SeqIO]]&,
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


cdef class Hashtable:
    cdef shared_ptr[CpHashtable] _ht_this

    cpdef bytes sanitize_kmer(self, object kmer)
    cpdef bytes sanitize_seq_kmer(self, object kmer)
    cdef HashIntoType sanitize_hash_kmer(self, object kmer) except -1
    cdef bytes _valid_sequence(self, str sequence)
    cdef CpKmer _build_kmer(self, object kmer) except *
    cdef FastxParserPtr _get_parser(self, object parser_or_filename) except *
    cdef list _get_raw_tables(self, uint8_t **, vector[uint64_t])


cdef class QFCounttable(Hashtable):
    cdef shared_ptr[CpQFCounttable] _qf_this


cdef class SmallCounttable(Hashtable):
    cdef shared_ptr[CpSmallCounttable] _st_this


cdef class Counttable(Hashtable):
    cdef shared_ptr[CpCounttable] _ct_this


cdef class CyclicCounttable(Hashtable):
    cdef shared_ptr[CpCyclicCounttable] _cct_this


cdef class Nodetable(Hashtable):
    cdef shared_ptr[CpNodetable] _nt_this


cdef class Hashgraph(Hashtable):
    cdef shared_ptr[CpHashgraph] _hg_this
    cdef SubsetPartition partitions
    # We keep an extra ref to the shared_ptr from partitions
    # to make sure dealloc ordering doesn't get clobbered
    cdef shared_ptr[CpSubsetPartition] partitions_ptr


cdef class Nodegraph(Hashgraph):
    cdef shared_ptr[CpNodegraph] _ng_this


cdef class Countgraph(Hashgraph):
    cdef shared_ptr[CpCountgraph] _cg_this


cdef class SmallCountgraph(Hashgraph):
    cdef shared_ptr[CpSmallCountgraph] _sg_this
