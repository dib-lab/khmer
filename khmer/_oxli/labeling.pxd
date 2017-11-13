from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.set cimport set
from libcpp.memory cimport unique_ptr, shared_ptr, weak_ptr
from libc.stdint cimport uint8_t, uint32_t, uint64_t, uintptr_t

from khmer._oxli.oxli_types cimport *
from khmer._oxli.hashing cimport Kmer, CpKmer, KmerSet, CpKmerFactory, CpKmerIterator
from khmer._oxli.graphs cimport (CpHashgraph, CpCountgraph, CpNodegraph, Hashgraph,
                                 Countgraph, Nodegraph)
from khmer._oxli.parsing cimport CpReadParser, CpSequence
from khmer._oxli.legacy_partitioning cimport (CpSubsetPartition, cp_pre_partition_info,
                                              SubsetPartition)
from khmer._oxli.utils cimport oxli_raise_py_error

cdef extern from "oxli/labelhash.hh" nogil:

    cdef cppclass CpLabelHash "oxli::LabelHash":
        CpLabeHash(CpHashgraph *)

        CpHashgraph * graph

        LabelSet all_labels

        size_t n_labels()
        void get_tag_labels(HashIntoType, LabelSet&) const
        void get_tags_from_label(Label, TagSet&) const
        void link_tag_and_label(const HashIntoType, const Label)
        unsigned int sweep_label_neighborhood(const string&,
                                              LabelSet&,
                                              unsigned int,
                                              bool, bool)
        void traverse_labels_and_resolve(const HashIntoTypeSet,
                                         LabelSet&)
        void save_labels_and_tags(string) except +oxli_raise_py_error
        void load_labels_and_tags(string) except +oxli_raise_py_error

        void label_across_high_degree_nodes(const char *,
                                            HashIntoTypeSet&,
                                            const Label)
        void get_labels_for_sequence(string&, LabelSet&) const

        void consume_seqfile_and_tag_with_labels[SeqIO](string &,
                                                       unsigned int &,
                                                       unsigned long long &) except +oxli_raise_py_error

        void consume_seqfile_and_tag_with_labels_readparser "consume_seqfile_and_tag_with_labels" [SeqIO](shared_ptr[CpReadParser[SeqIO]]&,
                                                        unsigned int &,
                                                        unsigned long long &) except +oxli_raise_py_error
        void consume_partitioned_fasta_and_tag_with_labels[SeqIO](string &,
                                                       unsigned int &,
                                                       unsigned long long &) except +oxli_raise_py_error
        void consume_sequence_and_tag_with_labels(const string &,   
                                                  unsigned long long&,
                                                  Label)
        void consume_sequence_and_tag_with_labels(const string &,   
                                                  unsigned long long&,
                                                  Label,
                                                  HashIntoTypeSet *)

cdef class GraphLabels:
    cdef shared_ptr[CpLabelHash] _lh_this
    cdef public Hashgraph graph
