from cython.operator cimport dereference as deref
from libcpp.memory cimport make_shared, shared_ptr

from khmer._oxli.graphs cimport Nodegraph, Hashgraph
from khmer._oxli.hashset cimport HashSet
from khmer._oxli.utils cimport _bstring
from khmer._oxli.utils import get_n_primes_near_x
from khmer._oxli.parsing cimport CpFastxReader

cdef class GraphLabels:

    def __cinit__(self, Hashgraph graph not None, *args, **kwargs):
        self.graph = graph
        self._lh_this = make_shared[CpLabelHash](graph._hg_this.get())

    @property
    def _default_sweep_radius(self):
        return (2 * self.graph.tag_density) + 1

    def consume_seqfile_and_tag_with_labels(self, str filename):
        cdef unsigned long long n_consumed = 0
        cdef unsigned int total_reads = 0
        cdef string _filename = _bstring(filename)
        deref(self._lh_this).consume_seqfile_and_tag_with_labels[CpFastxReader](_filename,
                                                                             total_reads,
                                                                             n_consumed)
        return total_reads, n_consumed

    def sweep_label_neighborhood(self, str sequence, object radius=0,
                                       bool break_on_stoptags=False,
                                       bool stop_big_traversals=False):
        cdef string _sequence = self.graph._valid_sequence(sequence)
        cdef unsigned int _radius = 0
        if radius is None:
            _radius = self._default_sweep_radius
        else:
            _radius = <unsigned int>radius
        cdef LabelSet found_labels
        deref(self._lh_this).sweep_label_neighborhood(_sequence, found_labels,
                                                      _radius, break_on_stoptags,
                                                      stop_big_traversals)
        
        cdef Label label
        for label in found_labels:
            yield label

    def consume_partitioned_fasta_and_tag_with_labels(self, str filename):
        cdef string _filename = _bstring(filename)
        cdef unsigned long long n_consumed = 0
        cdef unsigned int total_reads = 0
        deref(self._lh_this).consume_partitioned_fasta_and_tag_with_labels[CpFastxReader](_filename,
                                                                                       total_reads,
                                                                                       n_consumed)
        return total_reads, n_consumed

    def sweep_tag_neighborhood(self, str sequence, object radius=0,
                                       bool break_on_stoptags=False,
                                       bool stop_big_traversals=False):
        cdef string _sequence = self.graph._valid_sequence(sequence)
        cdef unsigned int _radius = 0
        if radius is None:
            _radius = self._default_sweep_radius
        else:
            _radius = <unsigned int>radius

        cdef HashSet tagged_kmers = HashSet(self.graph.ksize())
        deref(deref(self.graph._hg_this).partition).sweep_for_tags(_sequence,
                                                                   tagged_kmers.hs,
                                                                   deref(self.graph._hg_this).all_tags,
                                                                   _radius,
                                                                   break_on_stoptags,
                                                                   stop_big_traversals)
        return tagged_kmers

    def get_tag_labels(self, object tag):
        cdef HashIntoType _tag = self.graph.sanitize_hash_kmer(tag)
        cdef LabelSet labels
        deref(self._lh_this).get_tag_labels(_tag, labels)
        cdef Label label
        for label in labels:
            yield label

    def get_labels_for_sequence(self, str sequence):
        '''Return labels for all the tags found in this sequence.'''
        cdef string _sequence = self.graph._valid_sequence(sequence)
        cdef LabelSet labels
        deref(self._lh_this).get_labels_for_sequence(_sequence, labels)
        cdef Label label
        for label in labels:
            yield label

    def consume_sequence_and_tag_with_labels(self, str sequence, unsigned long long label):
        cdef string _sequence = _bstring(sequence)
        cdef unsigned long long n_consumed = 0
        deref(self._lh_this).consume_sequence_and_tag_with_labels(_sequence,
                                                               n_consumed, 
                                                               label)
        return n_consumed

    @property
    def n_labels(self):
        return deref(self._lh_this).n_labels()

    def labels(self):
        cdef Label label
        for label in deref(self._lh_this).all_labels:
            yield label

    def tags(self):
        '''Get all tagged k-mers as DNA strings.'''
        cdef HashIntoType st
        for st in deref(self.graph._hg_this).all_tags:
            yield deref(self.graph._hg_this).unhash_dna(st)

    def add_tag(self, object kmer):
        cdef HashIntoType _kmer = self.graph.sanitize_hash_kmer(kmer)
        deref(self.graph._hg_this).add_tag(_kmer)

    def link_tag_and_label(self, object tag, Label label):
        cdef HashIntoType _tag = self.graph.sanitize_hash_kmer(tag)
        deref(self._lh_this).link_tag_and_label(_tag, label)

    def label_across_high_degree_nodes(self, str sequence, HashSet hdns,
                                             unsigned long long label=0):
        cdef string _sequence = self.graph._valid_sequence(sequence)
        deref(self._lh_this).label_across_high_degree_nodes(_sequence.c_str(), hdns.hs, label)

    def assemble_labeled_path(self, *args, **kwargs):
        raise NotImplementedError()

    def save_labels_and_tags(self, str filename):
        deref(self._lh_this).save_labels_and_tags(_bstring(filename))

    def load_labels_and_tags(self, str filename):
        deref(self._lh_this).load_labels_and_tags(_bstring(filename))

    @staticmethod
    def load(str filename, Hashgraph graph):
        cdef GraphLabels gl = GraphLabels(graph)
        deref(gl._lh_this).load_labels_and_tags(_bstring(filename))
        return gl

    @staticmethod
    def NodeGraphLabels(int k, uint64_t starting_size, int n_tables, primes=[]):
        cdef vector[uint64_t] _primes
        if primes:
            _primes = primes
        else:
            _primes = get_n_primes_near_x(n_tables, starting_size)
        
        cdef Nodegraph graph = Nodegraph(k, starting_size, n_tables, primes=primes)
        return GraphLabels(graph)

    @staticmethod
    def CountGraphLabels(int k, uint64_t starting_size, int n_tables, primes=[]):
        cdef vector[uint64_t] _primes
        if primes:
            _primes = primes
        else:
            _primes = get_n_primes_near_x(n_tables, starting_size)
        
        cdef Countgraph graph = Countgraph(k, starting_size, n_tables, primes=primes)
        return GraphLabels(graph)

