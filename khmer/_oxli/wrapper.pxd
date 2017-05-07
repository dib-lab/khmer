from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.queue cimport queue
from libcpp.memory cimport unique_ptr, weak_ptr, shared_ptr
from libcpp.utility cimport pair
from libc.stdint cimport uint32_t, uint8_t, uint16_t, uint64_t, uintptr_t 


########################################################################
#
# Core: typedefs from oxli.hh.
#
########################################################################

cdef extern from "oxli/oxli.hh" namespace "oxli":
    ctypedef unsigned long long int HashIntoType
    ctypedef unsigned long long int Label
    ctypedef set[Label] LabelSet
    ctypedef set[HashIntoType] TagSet
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


cdef extern from "oxli/alphabets.hh" namespace "oxli":
    cdef string DNA_SIMPLE "oxli::alphabets::DNA_SIMPLE"
    cdef string DNAN_SIMPLE "oxli::alphabets::DNAN_SIMPLE"
    cdef string RNA_SIMPLE "oxli::alphabets::RNA_SIMPLE"
    cdef string RNAN_SIMPLE "oxli::alphabets::RNAN_SIMPLE"
    cdef string IUPAC_NUCL "oxli::alphabets::IUPAC_NUCL"
    cdef string IUPAC_AA "oxli::alphabets::IUPAC_AA"


########################################################################
#
# ReadParser: read parsing stuff, Read object
#
########################################################################

# C++ ostream wrapper code stolen shamelessly from stackoverflow
# http://stackoverflow.com/questions/30984078/cython-working-with-c-streams
# We need ostream to wrap ReadParser

cdef extern from "<iostream>" namespace "std":
    cdef cppclass ostream:
        ostream& write(const char*, int) except +

# obviously std::ios_base isn't a namespace, but this lets
# Cython generate the connect C++ code
cdef extern from "<iostream>" namespace "std::ios_base":
    cdef cppclass open_mode:
        pass
    cdef open_mode binary
    # you can define other constants as needed

cdef extern from "<fstream>" namespace "std":
    cdef cppclass ofstream(ostream):
        # constructors
        ofstream(const char*) except +
        ofstream(const char*, open_mode) except+

cdef extern from  "oxli/read_parsers.hh" namespace "oxli::read_parsers":
    cdef cppclass CpSequence "oxli::read_parsers::Read":
        string name
        string description
        string sequence
        string quality
        string cleaned_seq

        void reset()
        void write_fastx(ostream&)
        void set_cleaned_seq()        

    ctypedef pair[CpSequence,CpSequence] CpSequencePair \
        "oxli::read_parsers::ReadPair"

    cdef cppclass CpReadParser "oxli::read_parsers::ReadParser" [SeqIO]:
        CpReadParser(unique_ptr[SeqIO])
        CpReadParser(CpReadParser&)

        CpSequence get_next_read()
        CpSequencePair get_next_read_pair()
        CpSequencePair get_next_read_pair(uint8_t)

        uintptr_t get_num_reads()
        bool is_complete()
        void close()

    cdef cppclass CpFastxReader "oxli::read_parsers::FastxReader":
        CpFastxReader()
        CpFastxReader(const string&)
        CpFastxReader(CpFastxReader&)

        CpSequence get_next_read()
        bool is_complete()
        uintptr_t get_num_reads()
        void close()

    unique_ptr[CpReadParser[SeqIO]] get_parser[SeqIO](const string&) 
    ctypedef unique_ptr[CpReadParser[CpFastxReader]] FastxParserPtr

#
# Hashtable: Bindings for the existing CPython Hashtable wrapper.
#

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


########################################################################
#
# Traversal: wrapper for traversal.hh.
#
########################################################################


cdef extern from "oxli/traversal.hh" namespace "oxli":
    cdef cppclass CpTraverser "oxli::Traverser":
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


cdef extern from "oxli/assembler.hh" namespace "oxli":
    cdef cppclass CpLinearAssembler "oxli::LinearAssembler":
        CpLinearAssembler(CpHashgraph *)
    
        string assemble(const CpKmer, const CpHashgraph *) const
        string assemble_left(const CpKmer, const CpHashgraph *) const     
        string assemble_right(const CpKmer, const CpHashgraph *) const

        string assemble(const CpKmer) const
        string assemble_left(const CpKmer) const     
        string assemble_right(const CpKmer) const

    cdef cppclass CpCompactingAssembler(CpLinearAssembler):
        CpCompactingAssembler(CpHashgraph *)

    cdef cppclass CpSimpleLabeledAsssembler "oxli::SimpleLabeledAsssembler":
        CpSimpleLabeledAssembler(const CpLabelHash *)

        vector[string] assemble(const CpKmer)
        vector[string] assemble(const CpKmer, const CpHashgraph *) const

    cdef cppclass CpJunctionCountAssembler "oxli::JunctionCountAssembler":
        CpJunctionCountAssembler(CpHashgraph *)

        vector[string] assemble(const CpKmer) const
        vector[string] assemble(const CpKmer, const CpHashtable *) const
        uint16_t consume(string)
        void count_junction(CpKmer, CpKmer)
        BoundedCounterType get_junction_count(CpKmer, CpKmer) const


cdef extern from "oxli/hllcounter.hh" namespace "oxli":
    cdef cppclass CpHLLCounter "oxli::HLLCounter":
        CpHLLCounter(double, WordLength)
        CpHLLCounter(int, WordLength)

        void add(const string &)
        unsigned int consume_string(const string &)
        void consume_seqfile[SeqIO](const string &,
                                    bool,
                                    unsigned int &,
                                    unsigned long long &)

        void consume_seqfile[SeqIO](unique_ptr[CpReadParser[SeqIO]]&,
                                    bool,
                                    unsigned int &,
                                    unsigned long long &)
        unsigned int check_and_process_read(string &, bool &)
        bool check_and_normalize_read(string &) const
        uint64_t estimate_cardinality()
        void merge(CpHLLCounter &)
        double get_alpha()
        int get_p()
        int get_m()
        void set_ksize(WordLegth)
        int get_ksize()
        vector[int] get_M()
        double get_erate()
        void set_erate(double)


cdef extern from "oxli/partitioning.hh" namespace "oxli":

    ctypedef vector[HashIntoType] TagVector

    cdef cppclass CpComponent "oxli::Component":
        CpComponent()
        CpComponent(uint64_t)

        const uint64_t component_id
        vector[HashIntoType] tags

        void kill()
        bool is_alive() const

        void add_tag(HashIntoType)
        void add_tags(TagVector&)

        uint64_t get_n_tags() const
        uint64_t get_n_created() const
        uint64_t get_n_destroyed() const

    ctypedef shared_ptr[CpComponent] ComponentPtr
    ctypedef set[ComponentPtr] ComponentPtrSet
    ctypedef vector[ComponentPtr] ComponentPtrVector

    cdef cppclass CpGuardedHashCompMap "oxli::GuardedHashCompMap":
        map[HashIntoType, ComponentPtr] data

        ComponentPtr get(HashIntoType)
        void set(HashIntoType, ComponentPtr)
        bool contains(HashIntoType)

    cdef cppclass CpComponentMap "oxli::ComponentMap":
        CpComponentMap(WordLength, WordLength, uint64_t)

        void create_component(TagVector&)
        uint32_t create_and_merge_components(TagVector&)
        void map_tags_to_component(TagVector&, ComponentPtr&)
        uint32_t merge_components(ComponentPtr&, ComponentPtrSet&)

        bool contains(HashIntoType)
        ComponentPtr get(HashIntoType) const

        uint64_t get_n_components() const
        uint64_t get_n_tags() const
        weak_ptr[ComponentPtrVector] get_components()
        weak_ptr[CpGuardedHashCompMap] get_tag_component_map()

    cdef cppclass CpStreamingPartitioner "oxli::StreamingPartitioner":
        CpStreamingPartitioner(CpHashgraph * ) except +MemoryError
        CpStreamingPartitioner(CpHashgraph *, uint32_t) except +MemoryError
        
        CpHashgraph * graph
        uint64_t consume(string&) nogil except +MemoryError
        uint64_t  consume_pair(string&, string&) nogil except +MemoryError
        uint64_t consume_fasta(string&) except +MemoryError

        uint64_t seed_sequence(string&, TagVector&, KmerQueue&,
                           set[HashIntoType]&) except +MemoryError

        void find_connected_tags(queue[CpKmer]&, 
                                 TagVector&,
                                 set[HashIntoType]&) except +MemoryError

        void find_connected_tags(queue[CpKmer]&, 
                                 TagVector&,
                                 set[HashIntoType]&,
                                 bool) except +MemoryError


        uint64_t get_n_consumed() const
        uint32_t get_tag_density() const
        uint64_t get_n_components() const
        uint64_t get_n_tags() const

        ComponentPtr get_tag_component(string&) const
        ComponentPtr get_nearest_component(string&) const
        ComponentPtr get_nearest_component(CpKmer) const

        weak_ptr[ComponentPtrVector] get_components()
        weak_ptr[CpGuardedHashCompMap] get_tag_component_map()

