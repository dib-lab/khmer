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
    ctypedef unsigned long long int Label
    ctypedef set[Label] LabelSet
    ctypedef set[HashIntoType TagSet
    ctypedef set[HashIntoType] HashIntoTypeSet
    ctypedef unsigned char WordLength
    ctypedef unsigned short int BoundedCounterType
    ctypedef queue[CpKmer] KmerQueue
    ctypedef set[CpKmer] KmerSet
    ctypedef bool (*KmerFilter) (CpKmer kmer)
    ctypedef void (*CallbackFn)(const char *, void *, uint64_t, uint64_t)

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
    string _revcomp(const string&)


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

    cdef cppclass CpIParser "khmer::read_parsers::IParser":
        CpIParser()
        static CpIParser * const get_parser(string const &)
        bool is_complete()
        CpSequence get_next_read()
        CpSequencePair get_next_read_pair(uint8_t)
        size_t get_num_reads()

    cdef cppclass CpFastxParser(CpIParser) "khmer::read_parsers::FastxParser":
        CpFastxParser(const char *)


########################################################################
#
# Hashtable: Bindings for the existing CPython Hashtable wrapper.
#
########################################################################

# All we really need are the PyObject struct definitions
# for our extension objects.
cdef extern from "_khmer.hh":

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


cdef extern from "hashtable.hh" namespace "khmer":
    cdef cppclass CpHashtable "khmer::Hashtable":
        const WordLength ksize() const
        HashIntoType hash_dna(const char *) const
        HashIntoType hash_dna_top_strand(const char *) const
        HashIntpType hash_dna_bottom_strand(const char *) const
        string unhash_dna(HashIntoType) const
        void count(const char *)
        void count(HashIntoType)
        const BoundedCounterType get_count(const char *) const
        const BoundedCounterType get_count(HashIntoType) const
        void save(string)
        void load(string)
        uint32_t consume_string(const string &)
        bool check_and_normalize_read(string &) const
        uint32_t check_and_process_read(string &, bool &)
        void consume_fasta(string const &, uint32_t &, uint64_t &)
        void consume_fasta(CpIParser *, uint32_t &, uint64_t &)
        void set_use_bigcount(bool)
        bool get_use_bigcount()
        bool median_at_least(const string &, uint32_t cutoff)
        void get_median_count(const string &, BoundedCounterType &,
                              float &, float &)
        const uint64_t n_unique_kmers() const
        const uint64_t n_occupied() const
        vector[uint64_t] get_tablesizes() const
        const size_t n_tables() const
        void get_kmers(const string &, vector[string] &)
        void get_kmer_hashes(const string &, vector[HashIntoType] &) const
        void get_kmer_hashes_as_hashset(const string &, 
                                        set[HashIntoType]) const
        void get_kmer_counts(const string &, 
                             vector[BoundedCounterType] &) const
        BoundedCounterType get_min_count(const string &)
        BoundedCounterType get_max_count(const string &)
        uint64_t * abundance_distribution(CpIParser *, CpHashtable *)
        uint64_t * abundance_distribution(string, CpHashtable *)
        uint64_t trim_on_abundance(string, BoundedCounterType) const
        uint64_t trim_below_abundance(string, BoundedCounterType) const
        vector[uint32_t] find_spectral_error_positions(string, 
                                                       BoundedCounterType)

    cdef cppclass CpCounttable(CpHashtable) "khmer::Counttable":
        CpCounttable(WordLength, vector[uint64_t])

    cdef cppclass CpNodetable(CpHashtable) "khmer::Nodetable":
        CpNodetable(WordLength, vector[uint64_t])

cdef extern from "hashgraph.hh" namespace "khmer":
    cdef cppclass CpHashgraph(CpHashtable) "khmer::Hashgraph":
        uint32_t traverse_from_kmer(CpKmer, uint32_t, KmerSet&, uint32_t)
        void extract_unique_paths(string, uint32_t, float, vector[string])
        void calc_connected_graph_size(CpKmer, uint64_t&, KmerSet&,
                                       const uint64_t, bool) const
        uint32_t kmer_degree(HashIntoType, HashIntoType)
        uint32_t kmer_degree(const char *)
        void find_high_degree_nodes(const char *, set[HashIntoType] &) const

    cdef cppclass CpCountgraph(CpHashgraph) "khmer::Countgraph":
        CpCountgraph(WordLength, vector[uint64_t])

    cdef cppclass CpNodegraph(CpNodegraph) "khmer::Nodegraph":
        CpNodegraph(WordLength, vector[uint64_t])

cdef extern from "labelhash.hh" namespace "khmer":
    cdef cppclass CpLabelHash "khmer::LabelHash":
        CpLabelHash(CpHashgraph *)
        size_t n_labels() const
        void consume_fasta_and_tag_with_labels(string const &,
                                               uint32_t &,
                                               uint64_t &,
                                               CallbackFn,
                                               void *)
         void consume_fasta_and_tag_with_labels(string const &,
                                               uint32_t &,
                                               uint64_t &)
        void consume_fasta_and_tag_with_labels(CpIParser *,
                                               uint32_t &,
                                               uint64_t &,
                                               CallbackFn,
                                               void *)
         void consume_fasta_and_tag_with_labels(CpIParser *,
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

cdef CpHashtable * get_hashtable_ptr(object graph)


########################################################################
#
# Traversal: wrapper for traversal.hh.
#
########################################################################


cdef extern from "traversal.hh":
    cdef cppclass CpTraverser "khmer::Traverser":
        CpTraverser(CpHashgraph *)

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
# Assembler: wrapper for assembler.hh.
#
########################################################################


cdef extern from "assembler.hh" namespace "khmer":
    cdef cppclass CpLinearAssembler "khmer::LinearAssembler":
        CpLinearAssembler(CpHashgraph *)
    
        string assemble(const CpKmer, const CpHashgraph *) const
        string assemble_left(const CpKmer, const CpHashgraph *) const     
        string assemble_right(const CpKmer, const CpHashgraph *) const

        string assemble(const CpKmer) const
        string assemble_left(const CpKmer) const     
        string assemble_right(const CpKmer) const

