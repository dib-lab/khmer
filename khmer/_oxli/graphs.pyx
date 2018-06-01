from math import log

from cython.operator cimport dereference as deref
from cpython.buffer cimport (PyBuffer_FillInfo, PyBUF_FULL_RO)
from libc.stdint cimport uint64_t
from libc.stdint cimport uintptr_t as size_t

from libcpp.memory cimport shared_ptr, make_shared
from libcpp.vector cimport vector
from libcpp.set cimport set
from libcpp.string cimport string

from khmer._oxli.utils cimport _bstring, is_str, is_num
from khmer._oxli.utils import get_n_primes_near_x
from khmer._oxli.parsing cimport (CpFastxReader, CPyReadParser_Object,
                                  get_parser, CpReadParser, FastxParser,
                                  FastxParserPtr)
from khmer._oxli.hashset cimport HashSet
from khmer._oxli.legacy_partitioning cimport (CpSubsetPartition, SubsetPartition,
                                   cp_pre_partition_info, PrePartitionInfo)
from khmer._oxli.oxli_types cimport MAX_BIGCOUNT, HashIntoType
from khmer._oxli.traversal cimport Traverser

from khmer._khmer import ReadParser

CYTHON_TABLES = (Hashtable, Nodetable, Counttable, CyclicCounttable,
                 SmallCounttable,
                 QFCounttable, Nodegraph, Countgraph, SmallCountgraph)


cdef class Hashtable:

    cpdef bytes sanitize_seq_kmer(self, object kmer):
        '''Length sanitize a string k-mer and return as bytes.'''
        if len(kmer) != self.ksize():
            raise ValueError("Expected k-mer length {}"
                             " but got {}.".format(self.ksize(), len(kmer)))
        return _bstring(kmer)

    cpdef bytes sanitize_kmer(self, object kmer):
        '''Type and length sanitize a k-mer and return as bytes, reverse
        hashing if necessary.'''
        cdef bytes handled
        if is_num(kmer):
            handled = deref(self._ht_this).unhash_dna(<HashIntoType>kmer)
        elif isinstance(kmer, Kmer):
            handled = _bstring(kmer.kmer)
        elif is_str(kmer):
            handled = self.sanitize_seq_kmer(kmer)
        else:
            self._kmer_type_error(kmer)
        return handled

    cdef HashIntoType sanitize_hash_kmer(self, object kmer) except -1:
        '''Sanitize a hashed k-mer, or hash if not already.'''
        cdef HashIntoType handled
        if is_num(kmer):
            handled = <HashIntoType>kmer
        elif isinstance(kmer, Kmer):
            handled = kmer.kmer_u
        else:
            handled = self.hash(kmer)
        return handled

    cdef bytes _valid_sequence(self, str sequence):
        """Validate sequence argument and convert it to bytes"""
        if len(sequence) < self.ksize():
            raise ValueError("sequence length ({}) must >= the hashtable "
                             "k-mer size ({})".format(len(sequence),
                                                      self.ksize()))
        return _bstring(sequence)

    cdef CpKmer _build_kmer(self, object kmer) except *:
        '''Build a liboxli Kmer (CpKmer) from a hash or string.'''
        if type(kmer) is Kmer:
            return deref((<Kmer>kmer)._this.get())

        cdef bytes temp = self.sanitize_kmer(kmer)
        return deref(self._ht_this).build_kmer(temp)

    def _kmer_type_error(self, object kmer):
        raise TypeError("Object of type {0} can not be interpretted as "
                        " a k-mer".format(type(kmer)))

    def count(self, object kmer):
        """Increment the count of this k-mer.

        Synonym for 'add'.
        """
        self.add(kmer)

    def add(self, object kmer):
        """Increment the count of this k-mer

        `kmer` can be either a string or an integer representing the hashed
        value of the kmer.
        """
        if is_str(kmer):
            return deref(self._ht_this).add(self.sanitize_seq_kmer(kmer))
        elif is_num(kmer):
            return deref(self._ht_this).add(<HashIntoType>kmer)
        else:
            self._kmer_type_error(kmer)

    def hash(self, object kmer):
        """Compute the hash of this k-mer."""
        if is_num(kmer):
            return kmer

        data = self.sanitize_seq_kmer(kmer)
        return deref(self._ht_this).hash_dna(data)

    def reverse_hash(self, HashIntoType kmer_hash):
        """Turn a k-mer hash back into a DNA k-mer, if possible."""
        return deref(self._ht_this).unhash_dna(kmer_hash)

    def get(self, object kmer):
        """Retrieve the count for the given k-mer.

        `kmer` can be either a string or an integer representing the hashed
        value of the kmer.

        For Nodetables and Counttables, this function will fail if the
        supplied k-mer contains non-ACGT characters.
        """
        if is_str(kmer):
            _kmer = self.sanitize_seq_kmer(kmer)
            return deref(self._ht_this).get_count(_kmer)
        elif is_num(kmer):
            return deref(self._ht_this).get_count(<HashIntoType> kmer)
        else:
            self._kmer_type_error(kmer)


    def ksize(self):
        """k-mer size"""
        return deref(self._ht_this).ksize()

    def hashsizes(self):
        """Size of hash tables used."""
        return deref(self._ht_this).get_tablesizes()

    def get_kmers(self, str sequence):
        """Generate an ordered list of all k-mers in sequence."""
        cdef bytes data = self._valid_sequence(sequence)
        cdef vector[string] kmers
        deref(self._ht_this).get_kmers(data, kmers)
        return kmers

    def consume(self, str sequence):
        """Increment the counts of all of the k-mers in the sequence."""
        cdef bytes data = self._valid_sequence(sequence)
        return deref(self._ht_this).consume_string(data)

    def get_kmer_counts(self, str sequence):
        """Retrieve an ordered list of the counts of all k-mers in sequence."""
        cdef bytes data = self._valid_sequence(sequence)
        cdef vector[BoundedCounterType] counts
        deref(self._ht_this).get_kmer_counts(data, counts)
        return counts

    def get_min_count(self, str sequence):
        """Get the smallest count of all the k-mers in the string."""
        cdef bytes data = self._valid_sequence(sequence)
        return deref(self._ht_this).get_min_count(data)

    def get_max_count(self, str sequence):
        """Get the larget count of all the k-mers in the string."""
        cdef bytes data = self._valid_sequence(sequence)
        return deref(self._ht_this).get_max_count(data)

    def get_median_count(self, str sequence):
        """median, average, and stddev of the k-mer counts in sequence."""
        cdef bytes data = self._valid_sequence(sequence)
        cdef BoundedCounterType med = 0
        cdef float average = 0
        cdef float stddev = 0

        deref(self._ht_this).get_median_count(data, med, average, stddev)
        return (med, average, stddev)

    def median_at_least(self, str sequence, int median):
        '''Check if median k-mer count is at least the given value.'''
        cdef bytes data = self._valid_sequence(sequence)
        return <bool>deref(self._ht_this).median_at_least(data, median)

    def get_kmer_hashes(self, str sequence):
        """Retrieve hashes of all k-mers in sequence.

        Hashes are returned in the same order as k-mers appear in sequence.
        """
        cdef bytes data = self._valid_sequence(sequence)
        cdef vector[HashIntoType] hashes
        deref(self._ht_this).get_kmer_hashes(data, hashes)
        return hashes

    def trim_on_abundance(self, str sequence, int abundance):
        """Trim sequence at first k-mer below the given abundance."""
        cdef bytes data = self._valid_sequence(sequence)
        trimmed_at = deref(self._ht_this).trim_on_abundance(data, abundance)
        return sequence[:trimmed_at], trimmed_at

    def trim_below_abundance(self, str sequence, int abundance):
        """Trim sequence at first k-mer above the given abundance."""
        cdef bytes data = self._valid_sequence(sequence)
        cdef int trimmed_at = deref(self._ht_this).trim_below_abundance(data, abundance)
        return sequence[:trimmed_at], trimmed_at

    def find_spectral_error_positions(self, str sequence, int max_count):
        """Identify positions of low-abundance k-mers."""
        cdef bytes data = self._valid_sequence(sequence)
        posns = (deref(self._ht_this).find_spectral_error_positions(data,
                                                                   max_count))
        return posns

    cdef FastxParserPtr _get_parser(self, object parser_or_filename) except *:
        cdef FastxParserPtr _parser
        if isinstance(parser_or_filename, FastxParser):
            _parser = (<FastxParser>parser_or_filename)._this
        elif isinstance(parser_or_filename, ReadParser):
            _parser = (<CPyReadParser_Object*>parser_or_filename).parser
        elif is_str(parser_or_filename):
            _parser = get_parser[CpFastxReader](_bstring(parser_or_filename))
        else:
            raise TypeError('argument does not appear to be a parser or a '
                            'filename: {}'.format(parser_or_filename))
        return _parser


    def consume_seqfile(self, object parser_or_filename):
        """Count all k-mers from file_name."""
        cdef unsigned long long n_consumed = 0
        cdef unsigned int total_reads = 0
        cdef FastxParserPtr _parser = self._get_parser(parser_or_filename)
        with nogil:
            deref(self._ht_this).consume_seqfile[CpFastxReader](\
                _parser, total_reads, n_consumed
            )
        return total_reads, n_consumed

    def consume_seqfile_with_mask(self, object parser_or_filename,
                                  Hashtable mask, int threshold=0,
                                  bool consume_masked=False):
        cdef unsigned long long n_consumed = 0
        cdef unsigned int total_reads = 0
        cdef FastxParserPtr _parser = self._get_parser(parser_or_filename)
        cdef CpHashtable * cmask = mask._ht_this.get()
        with nogil:
            deref(self._ht_this).consume_seqfile_with_mask[CpFastxReader](\
                _parser, cmask, threshold, total_reads, n_consumed, consume_masked
            )
        return total_reads, n_consumed

    def consume_seqfile_banding(self, object parser_or_filename, int num_bands,
                                int band):
        """Count all k-mers from file_name."""
        cdef unsigned long long n_consumed = 0
        cdef unsigned int total_reads = 0
        cdef FastxParserPtr _parser = self._get_parser(parser_or_filename)
        with nogil:
            deref(self._ht_this).consume_seqfile_banding[CpFastxReader](\
                _parser, num_bands, band, total_reads, n_consumed
            )
        return total_reads, n_consumed

    def consume_seqfile_banding_with_mask(self, object parser_or_filename,
                                          int num_bands, int band,
                                          Hashtable mask, int threshold=0,
                                          bool consume_masked=False):
        cdef unsigned long long n_consumed = 0
        cdef unsigned int total_reads = 0
        cdef FastxParserPtr _parser = self._get_parser(parser_or_filename)
        cdef CpHashtable * cmask = mask._ht_this.get()
        with nogil:
            deref(self._ht_this).\
                consume_seqfile_banding_with_mask[CpFastxReader](\
                    _parser, num_bands, band, cmask, threshold, total_reads,
                    n_consumed, consume_masked
                )
        return total_reads, n_consumed

    def abundance_distribution(self, object parser_or_filename,
                               Hashtable tracking):
        """Calculate the k-mer abundance distribution over input reads."""
        cdef FastxParserPtr _parser = self._get_parser(parser_or_filename)
        cdef CpHashtable * _tracking = tracking._ht_this.get()
        cdef uint64_t * x
        with nogil:
            x = deref(self._ht_this).abundance_distribution[CpFastxReader](\
                _parser, _tracking
            )

        abunds = []
        for i in range(MAX_BIGCOUNT):
            abunds.append(x[i])
        return abunds

    def save(self, file_name):
        """Save the graph to the specified file."""
        deref(self._ht_this).save(_bstring(file_name))

    @classmethod
    def load(cls, file_name):
        """Load the graph from the specified file."""
        cdef Hashtable table = cls(1, 1, 1)
        deref(table._ht_this).load(_bstring(file_name))
        return table

    def n_unique_kmers(self):
        """Estimate of the number of unique kmers stored."""
        return deref(self._ht_this).n_unique_kmers()

    def n_occupied(self):
        """Estimate of the number of occupied slots in the storage."""
        return deref(self._ht_this).n_occupied()

    def n_tables(self):
        """Number of tables used in the storage."""
        return deref(self._ht_this).n_tables()

    def set_use_bigcount(self, bigcount):
        deref(self._ht_this).set_use_bigcount(<bool>bigcount)

    def get_use_bigcount(self):
        return deref(self._ht_this).get_use_bigcount()

    def get_kmer_hashes_as_hashset(self, str sequence):
        cdef HashSet hashes = HashSet(self.ksize())
        deref(self._ht_this).get_kmer_hashes_as_hashset(_bstring(sequence),
                                                        hashes.hs)
        return hashes

    cdef list _get_raw_tables(self, uint8_t ** table_ptrs, vector[uint64_t] sizes):
        cdef Py_buffer buf_info
        cdef object view
        cdef list views = []
        for table_idx in range(0, len(sizes)):
            PyBuffer_FillInfo(&buf_info, None, table_ptrs[table_idx],
                              sizes[table_idx], 0, PyBUF_FULL_RO)
            view = PyMemoryView_FromBuffer(&buf_info)
            views.append(view)
        return views

    def get_raw_tables(self):
        cdef uint8_t ** table_ptrs = deref(self._ht_this).get_raw_tables()
        cdef vector[uint64_t] sizes = deref(self._ht_this).get_tablesizes()
        return self._get_raw_tables(table_ptrs, sizes)


cdef class QFCounttable(Hashtable):
    """Count kmers using a counting quotient filter.

    The counting quotient filter (CQF) is an extension of the quotient filter
    that supports counting in addition to simple membership testing. A CQF has
    better cache locality compared to (Small)Counttable which increases
    performance.

    Each new k-mer uses one slot, and the number of slots used per k-mer
    increases the more often the same k-mer is entered into the CQF. As a result
    the CQF can be "full" and will stop accepting calls to `add` and `count`.

    Parameters
    ----------
    k : integer
        k-mer size

    size : integer
        Set the number of slots used by the counting quotient filter. This
        determines the amount of memory used and how many k-mers can be entered
        into the datastructure. Each slot uses roughly 1.3 bytes.
    """

    def __cinit__(self, int k, uint64_t size):
        # size has to be a power of two
        power_of_two = ((size & (size - 1) == 0) and
                        (size != 0))
        if not power_of_two:
            raise ValueError("size has to be a power of two, not"
                             " {}.".format(size))
        if type(self) is QFCounttable:
            self._qf_this = make_shared[CpQFCounttable](k, <uint64_t>log(size, 2))
            self._ht_this = <shared_ptr[CpHashtable]>self._qf_this


    @classmethod
    def load(cls, file_name):
        """Load the graph from the specified file."""
        cdef QFCounttable table = cls(1, 1)
        deref(table._qf_this).load(_bstring(file_name))
        return table

cdef class Counttable(Hashtable):

    def __cinit__(self, int k, uint64_t starting_size, int n_tables,
                  primes=None):
        if primes is None:
            primes = list()
        cdef vector[uint64_t] _primes
        if type(self) is Counttable:
            if primes:
                _primes = primes
            else:
                _primes = get_n_primes_near_x(n_tables, starting_size)
            self._ct_this = make_shared[CpCounttable](k, _primes)
            self._ht_this = <shared_ptr[CpHashtable]>self._ct_this


cdef class CyclicCounttable(Hashtable):

    def __cinit__(self, int k, uint64_t starting_size, int n_tables,
                  primes=None):
        if primes is None:
            primes = list()
        cdef vector[uint64_t] _primes
        if type(self) is CyclicCounttable:
            if primes:
                _primes = primes
            else:
                _primes = get_n_primes_near_x(n_tables, starting_size)
            self._cct_this = make_shared[CpCyclicCounttable](k, _primes)
            self._ht_this = <shared_ptr[CpHashtable]>self._cct_this


cdef class SmallCounttable(Hashtable):

    def __cinit__(self, int k, uint64_t starting_size, int n_tables,
                  primes=None):
        if primes is None:
            primes = list()
        cdef vector[uint64_t] _primes
        if type(self) is SmallCounttable:
            if primes:
                _primes = primes
            else:
                _primes = get_n_primes_near_x(n_tables, starting_size)
            self._st_this = make_shared[CpSmallCounttable](k, _primes)
            self._ht_this = <shared_ptr[CpHashtable]>self._st_this

    def get_raw_tables(self):
        cdef uint8_t ** table_ptrs = deref(self._st_this).get_raw_tables()
        cdef vector[uint64_t] sizes = deref(self._st_this).get_tablesizes()
        for i in range(len(sizes)):
            sizes[i] = (sizes[i] // 2) + 1
        return self._get_raw_tables(table_ptrs, sizes)


cdef class Nodetable(Hashtable):

    def __cinit__(self, int k, uint64_t starting_size, int n_tables,
                  primes=None):
        if primes is None:
            primes = list()
        cdef vector[uint64_t] _primes
        if type(self) is Nodetable:
            if primes:
                _primes = primes
            else:
                _primes = get_n_primes_near_x(n_tables, starting_size)
            self._nt_this = make_shared[CpNodetable](k, _primes)
            self._ht_this = <shared_ptr[CpHashtable]>self._nt_this


cdef class Hashgraph(Hashtable):

    def __cinit__(self, *args, **kwargs):
        self.partitions = None

    @property
    def partition(self):
        if self.partitions is None:
            self.partitions = SubsetPartition(self)
            self.partitions._this = deref(self._hg_this).partition
            self.partitions_ptr = self.partitions._this
        return self.partitions

    def neighbors(self, object kmer):
        '''Get a list of neighbor nodes for this k-mer.'''
        cdef Traverser traverser = Traverser(self)
        return list(traverser._neighbors(self._build_kmer(kmer)))

    def calc_connected_graph_size(self, object kmer, max_size=0,
                                  break_on_circumference=False):
        '''Find the number of nodes connected to this k-mer.'''
        cdef CpKmer _kmer = self._build_kmer(_bstring(kmer))
        cdef unsigned long long _size = 0
        cdef uint64_t _max_size = max_size
        cdef bool _break = break_on_circumference
        cdef KmerSet keeper
        cdef CpHashgraph * ptr = self._hg_this.get() # need tmp ref for nogil

        with nogil:
            deref(ptr).calc_connected_graph_size(_kmer, _size,
                                                 keeper, _max_size,
                                                 _break)
        return _size

    def kmer_degree(self, object kmer):
        '''Calculate the number of immediate neighbors this k-mer has
        the graph.'''
        cdef bytes _kmer = self.sanitize_kmer(kmer)
        return deref(self._hg_this).kmer_degree(_kmer)

    def count_kmers_within_radius(self, object kmer, int radius, int max_count=0):
        '''Calculate the number of neighbors with given radius in the graph.'''
        cdef unsigned int n
        cdef uint32_t _radius = radius
        cdef uint32_t _max_count = max_count
        cdef CpKmer _kmer = self._build_kmer(kmer)
        cdef KmerSet seen
        cdef CpHashgraph * ptr = self._hg_this.get()
        with nogil:
            n = deref(ptr).traverse_from_kmer(_kmer, _radius,
                                              seen, _max_count)
        return n

    def find_high_degree_nodes(self, str sequence):
        '''Examine the given sequence for degree > 2 nodes and add to
        list; used in graph contraction.'''
        cdef HashSet hdns = HashSet(self.ksize())
        _sequence = self._valid_sequence(sequence)
        deref(self._hg_this).find_high_degree_nodes(_sequence,
                                                    hdns.hs)
        return hdns


    def traverse_linear_path(self, object kmer, HashSet hdns,
                             Nodegraph stop_filter=None):
        '''Traverse the path through the graph starting with the given
        k-mer and avoiding high-degree nodes, finding (and returning)
        traversed k-mers and any encountered high-degree nodes.'''
        cdef HashSet adj = HashSet(self.ksize())
        cdef HashSet visited = HashSet(self.ksize())
        cdef CpKmer _kmer = self._build_kmer(kmer)
        cdef CpNodegraph * _stop_filter = stop_filter._ng_this.get()
        cdef int size = deref(self._hg_this).traverse_linear_path(_kmer,
                                                                 adj.hs,
                                                                 visited.hs,
                                                                 deref(_stop_filter),
                                                                 hdns.hs)
        return size, adj, visited

    def extract_unique_paths(self, str sequence, unsigned int min_length, float
                             min_unique_f):
        cdef vector[string] results
        deref(self._hg_this).extract_unique_paths(_bstring(sequence), min_length,
                                                  min_unique_f, results)
        return results

    def consume_and_tag(self, str sequence):
        '''Consume a sequence and tag it.'''
        cdef unsigned long long n_consumed = 0
        deref(self._hg_this).consume_sequence_and_tag(_bstring(sequence),
                                                     n_consumed)
        return n_consumed

    def get_tags_and_positions(self, str sequence):
        '''Retrieve tags and their positions in a sequence.'''
        cdef list result = []
        cdef int pos
        cdef WordLength K = deref(self._hg_this).ksize()
        cdef HashIntoType kmer
        for pos in range(0, len(sequence)-K+1):
            kmer = deref(self._hg_this).hash_dna(_bstring(sequence[pos:pos+K]))
            if deref(self._hg_this).has_tag(kmer):
                result.append((pos+1, kmer))
        return result

    def get_tags_for_sequence(self, str sequence):
        '''Get the tags present in a sequence.'''
        cdef string _sequence = self._valid_sequence(sequence)
        cdef HashSet hs = HashSet(self.ksize())
        deref(self._hg_this).get_tags_for_sequence(_sequence, hs.hs)
        return hs

    def find_all_tags_list(self, object kmer):
        '''Find all tags within range of the given k-mer, return as list'''
        cdef CpKmer _kmer = self._build_kmer(kmer)
        cdef HashSet result = HashSet(self.ksize())
        cdef set[HashIntoType] * tags = &(result.hs)
        cdef shared_ptr[CpHashgraph] this = self._hg_this

        with nogil:
            deref(deref(self._hg_this).partition).find_all_tags(_kmer, deref(tags),
                                                                deref(this).all_tags)

        return result

    def consume_seqfile_and_tag(self, object parser_or_filename):
        '''Consume all sequences in a FASTA/FASTQ file and tag the resulting
        graph.'''
        cdef unsigned long long n_consumed = 0
        cdef unsigned int total_reads = 0
        cdef FastxParserPtr _parser = self._get_parser(parser_or_filename)

        with nogil:
            deref(self._hg_this).\
                consume_seqfile_and_tag_readparser[CpFastxReader](_parser,
                                                                  total_reads,
                                                                  n_consumed)

        return total_reads, n_consumed

    def print_tagset(self, str filename):
        '''Print out all of the tags.'''
        deref(self._hg_this).print_tagset(_bstring(filename))

    def add_tag(self, object kmer):
        '''Add a k-mer to the tagset.'''
        cdef HashIntoType _kmer = self.sanitize_hash_kmer(kmer)
        deref(self._hg_this).add_tag(_kmer)

    def get_tagset(self):
        '''Get all tagged k-mers as DNA strings.'''
        cdef HashIntoType st
        cdef list all_tags = []
        for st in deref(self._hg_this).all_tags:
            all_tags.append(deref(self._hg_this).unhash_dna(st))
        return all_tags

    def tags(self):
        '''Get all tagged k-mers as DNA strings.'''
        cdef HashIntoType st
        for st in deref(self._hg_this).all_tags:
            yield deref(self._hg_this).unhash_dna(st)

    def load_tagset(self, str filename, clear_tags=True):
        '''Load tags from a file.'''
        deref(self._hg_this).load_tagset(_bstring(filename), clear_tags)

    def save_tagset(self, str filename):
        '''Save tags to a file.'''
        deref(self._hg_this).save_tagset(_bstring(filename))

    @property
    def n_tags(self):
        '''Return the count of all tags.'''
        return deref(self._hg_this).n_tags()

    def divide_tags_into_subsets(self, int subset_size=0):
        '''Divide tags equally up into subsets of given size.'''
        cdef set[HashIntoType] divvy
        deref(self._hg_this).divide_tags_into_subsets(subset_size, divvy)
        cdef HashSet hs = HashSet(self.ksize())
        hs.hs = divvy
        return hs

    @property
    def tag_density(self):
        '''Get the tagging density.'''
        return deref(self._hg_this)._get_tag_density()

    @tag_density.setter
    def tag_density(self, int density):
        '''Set the tagging density.'''
        deref(self._hg_this)._set_tag_density(density)

    def do_subset_partition(self, object start_kmer, object end_kmer,
                                  bool break_on_stoptags=False,
                                  bool stop_big_traversals=False):
        '''Partition the graph starting from a given subset of tags.'''

        cdef SubsetPartition subset = SubsetPartition(self)
        cdef CpSubsetPartition * subset_ptr = subset._this.get()
        cdef HashIntoType start = self.sanitize_hash_kmer(start_kmer)
        cdef HashIntoType end = self.sanitize_hash_kmer(end_kmer)
        cdef bool cbreak = break_on_stoptags
        cdef bool cstop = stop_big_traversals

        with nogil:
            deref(subset_ptr).do_partition(start, end, cbreak, cstop)

        return subset


    def find_all_tags(self, object kmer):
        '''Starting from the given k-mer, find all closely connected tags.'''
        cdef CpKmer _kmer = self._build_kmer(kmer)
        cdef PrePartitionInfo ppi = PrePartitionInfo.create(_kmer)

        with nogil:
            deref(deref(self._hg_this).partition).find_all_tags(_kmer,
                                                                deref(ppi._this).tagged_kmers,
                                                                deref(self._hg_this).all_tags)
            deref(self._hg_this).add_kmer_to_tags(_kmer.kmer_u)

        return ppi


    def assign_partition_id(self, PrePartitionInfo ppi):
        '''Assign a partition ID to a given tag.'''
        cdef cp_pre_partition_info * cppi = ppi._this.get()
        cdef PartitionID pi
        pi = deref(deref(self._hg_this).partition).assign_partition_id(deref(cppi).kmer,
                                                                       deref(cppi).tagged_kmers)
        return pi

    def output_partitions(self, str filename, str output, bool
                                output_unassigned=False):
        '''Write out sequences in given filename to another file, annotating '''
        '''with partition IDs.'''
        n_partitions = deref(deref(self._hg_this).partition).\
                            output_partitioned_file(_bstring(filename),
                                                    _bstring(output),
                                                    output_unassigned)
        return n_partitions

    def load_partitionmap(self, str filename):
        '''Load a partitionmap for the master subset.'''
        deref(deref(self._hg_this).partition).load_partitionmap(_bstring(filename))

    def save_partitionmap(self, str filename):
        '''Save a partitionmap for the master subset.'''
        deref(deref(self._hg_this).partition).save_partitionmap(_bstring(filename))

    def _validate_partitionmap(self):
        '''Run internal validation checks.'''
        deref(deref(self._hg_this).partition)._validate_pmap()

    def consume_partitioned_fasta(self, filename):
        '''Count all k-mers in a given file'''
        cdef unsigned long long n_consumed = 0
        cdef unsigned int total_reads = 0
        cdef string _filename = _bstring(filename)
        deref(self._hg_this).consume_partitioned_fasta[CpFastxReader](_filename,
                                                                     total_reads,
                                                                     n_consumed)
        return total_reads, n_consumed

    def merge_subset(self, SubsetPartition subset):
        '''Merge the given subset into this one.'''
        deref(deref(self._hg_this).partition).merge(subset._this.get())

    def merge_subset_from_disk(self, str filename):
        '''Merge the given subset (filename) into this one.'''
        deref(deref(self._hg_this).partition).merge_from_disk(_bstring(filename))

    def count_partitions(self):
        '''Count the number of partitions in the master partitionmap.'''
        return self.partition.count_partitions()

    def set_partition_id(self, object kmer, PartitionID pid):
        '''Set the partition ID for this tag.'''
        cdef string start = self.sanitize_kmer(kmer)
        deref(deref(self._hg_this).partition).set_partition_id(start, pid)

    def join_partitions(self, PartitionID p1, PartitionID p2):
        '''Join the partitions of these two tags.'''
        return deref(deref(self._hg_this).partition).join_partitions(p1, p2)

    def get_partition_id(self, object kmer):
        '''Get the partition ID of this tag.'''
        cdef string _kmer = self.sanitize_kmer(kmer)
        return deref(deref(self._hg_this).partition).get_partition_id(_kmer)

    def repartition_largest_partition(self, Countgraph counts not None,
                                            unsigned int distance,
                                            unsigned int threshold,
                                            unsigned int frequency,
                                            SubsetPartition subs=None):
        '''Repartition the largest partition (in the face of stop tags).'''

        cdef shared_ptr[CpSubsetPartition] subs_ptr
        if subs is None:
            subs_ptr = deref(self._hg_this).partition
        else:
            subs_ptr = subs._this

        cdef unsigned long next_largest
        next_largest = deref(subs_ptr).\
                repartition_largest_partition(distance,
                                              threshold,
                                              frequency,
                                              deref(counts._cg_this))
        return next_largest

    def load_stop_tags(self, object filename, clear_tags=False):
        '''Load the set of stop tags.'''
        deref(self._hg_this).load_stop_tags(_bstring(filename), clear_tags)

    def save_stop_tags(self, object filename):
        '''Save the set of stop tags.'''
        deref(self._hg_this).save_stop_tags(_bstring(filename))

    def print_stop_tags(self, filename):
        '''Print out the set of stop tags.'''
        deref(self._hg_this).print_stop_tags(_bstring(filename))

    def trim_on_stoptags(self, str sequence):
        '''Trim the reads on the given stop tags.'''
        cdef size_t trim_at
        cdef CpHashgraph * ptr = self._hg_this.get()
        cdef string cseq = _bstring(sequence)
        with nogil:
            trim_at = deref(ptr).trim_on_stoptags(cseq)
        return sequence[:trim_at], trim_at

    def add_stop_tag(self, object kmer):
        '''Add this k-mer as a stop tag.'''
        cdef HashIntoType _kmer = self.sanitize_hash_kmer(kmer)
        deref(self._hg_this).add_stop_tag(_kmer)

    def get_stop_tags(self):
        '''Return a DNA list of all of the stop tags.'''
        cdef HashIntoType st
        cdef list stop_tags = []
        for st in deref(self._hg_this).stop_tags:
            stop_tags.append(deref(self._hg_this).unhash_dna(st))
        return stop_tags

    def iter_stop_tags(self):
        '''Return a DNA list of all of the stop tags.'''
        cdef HashIntoType st
        for st in deref(self._hg_this).stop_tags:
            yield deref(self._hg_this).unhash_dna(st)


cdef class Countgraph(Hashgraph):

    def __cinit__(self, int k, uint64_t starting_size, int n_tables,
                  primes=None):
        if primes is None:
            primes = list()
        cdef vector[uint64_t] _primes
        if type(self) is Countgraph:
            if primes:
                _primes = primes
            else:
                _primes = get_n_primes_near_x(n_tables, starting_size)
            self._cg_this = make_shared[CpCountgraph](k, _primes)
            self._hg_this = <shared_ptr[CpHashgraph]>self._cg_this
            self._ht_this = <shared_ptr[CpHashtable]>self._hg_this

    def do_subset_partition_with_abundance(self, BoundedCounterType min_count,
                                                 BoundedCounterType max_count,
                                                 object start_kmer=0,
                                                 object end_kmer=0,
                                                 bool break_on_stop_tags=False,
                                                 bool stop_big_traversals=False):

        cdef HashIntoType _start_kmer = self.sanitize_hash_kmer(start_kmer)
        cdef HashIntoType _end_kmer = self.sanitize_hash_kmer(end_kmer)
        cdef bool _break_on_stop_tags = break_on_stop_tags
        cdef bool _stop_big_traversals = stop_big_traversals
        cdef SubsetPartition subset = SubsetPartition(self)
        cdef shared_ptr[CpSubsetPartition] subset_ptr = subset._this

        with nogil:
            deref(subset_ptr).do_partition_with_abundance(_start_kmer,
                                                          _end_kmer,
                                                          min_count,
                                                          max_count,
                                                          break_on_stop_tags,
                                                          stop_big_traversals)

        return subset


cdef class SmallCountgraph(Hashgraph):

    def __cinit__(self, int k, uint64_t starting_size, int n_tables,
                  primes=None):
        if primes is None:
            primes = list()
        cdef vector[uint64_t] _primes
        if type(self) is SmallCountgraph:
            if primes:
                _primes = primes
            else:
                _primes = get_n_primes_near_x(n_tables, starting_size)
            self._sg_this = make_shared[CpSmallCountgraph](k, _primes)
            self._hg_this = <shared_ptr[CpHashgraph]>self._sg_this
            self._ht_this = <shared_ptr[CpHashtable]>self._hg_this

    def get_raw_tables(self):
        cdef uint8_t ** table_ptrs = deref(self._sg_this).get_raw_tables()
        cdef vector[uint64_t] sizes = deref(self._sg_this).get_tablesizes()
        for i in range(len(sizes)):
            sizes[i] = sizes[i] // 2 + 1
        return self._get_raw_tables(table_ptrs, sizes)



cdef class Nodegraph(Hashgraph):

    def __cinit__(self, int k, uint64_t starting_size, int n_tables,
                  primes=None):
        if primes is None:
            primes = list()
        cdef vector[uint64_t] _primes
        if type(self) is Nodegraph:
            if primes:
                _primes = primes
            else:
                _primes = get_n_primes_near_x(n_tables, starting_size)
            self._ng_this = make_shared[CpNodegraph](k, _primes)
            self._hg_this = <shared_ptr[CpHashgraph]>self._ng_this
            self._ht_this = <shared_ptr[CpHashtable]>self._hg_this

    def update(self, Nodegraph other):
        deref(self._ng_this).update_from(deref(other._ng_this))
