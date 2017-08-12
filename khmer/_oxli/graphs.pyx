# cython: c_string_type=unicode, c_string_encoding=utf8
from math import log

from cython.operator cimport dereference as deref
from cpython.buffer cimport (PyBuffer_FillInfo, PyBUF_FULL_RO)
from libc.stdint cimport uint64_t
from libc.stdint cimport uintptr_t as size_t

from libcpp.memory cimport unique_ptr, shared_ptr, make_shared
from libcpp.vector cimport vector
from libcpp.set cimport set
from libcpp.string cimport string

from .utils cimport _bstring
from .utils import get_n_primes_near_x
from .parsing cimport (CpFastxReader, CPyReadParser_Object, get_parser,
                      CpReadParser, FastxParserPtr)
from .hashset cimport HashSet
from .oxli_types cimport MAX_BIGCOUNT, HashIntoType
from .traversal cimport Traverser

from .._khmer import Countgraph as PyCountgraph
from .._khmer import Nodegraph as PyNodegraph
from .._khmer import GraphLabels as PyGraphLabels
from .._khmer import ReadParser


CYTHON_TABLES = (Hashtable, Nodetable, Counttable, SmallCounttable,
                 QFCounttable, Nodegraph, Countgraph, SmallCountgraph)
CPYTHON_TABLES = (PyCountgraph, PyNodegraph)


cdef CpHashgraph * get_hashgraph_ptr(object graph):
    if not (isinstance(graph, PyCountgraph) or isinstance(graph, PyNodegraph)):
        return NULL

    cdef CPyHashgraph_Object* ptr = <CPyHashgraph_Object*> graph
    return deref(ptr).hashgraph


cdef CpLabelHash * get_labelhash_ptr(object labels):
    if not isinstance(labels, PyGraphLabels):
        return NULL

    cdef CPyGraphLabels_Object * ptr = <CPyGraphLabels_Object*> labels
    return deref(ptr).labelhash


cdef class Hashtable:
    def count(self, kmer):
        """Increment the count of this k-mer.

        Synonym for 'add'.
        """
        self.add(kmer)

    def add(self, kmer):
        """Increment the count of this k-mer

        `kmer` can be either a string or an integer representing the hashed
        value of the kmer.
        """
        if isinstance(kmer, basestring):
            temp = kmer.encode('utf-8')
            return deref(self._ht_this).add(<char*>temp)
        # assume kmer is an integer representing the hash value
        else:
            return deref(self._ht_this).add(<uint64_t>kmer)

    def hash(self, kmer):
        """Compute the hash of this k-mer."""
        if len(kmer) != self.ksize():
            raise ValueError("Expected k-mer length {}"
                             " but got {}.".format(self.ksize(), len(kmer)))
        data = _bstring(kmer)
        return deref(self._ht_this).hash_dna(data)

    def reverse_hash(self, kmer_hash):
        """Turn a k-mer hash back into a DNA k-mer, if possible."""
        return deref(self._ht_this).unhash_dna(kmer_hash)

    def get(self, kmer):
        """Retrieve the count for the given k-mer.

        `kmer` can be either a string or an integer representing the hashed
        value of the kmer.

        For Nodetables and Counttables, this function will fail if the
        supplied k-mer contains non-ACGT characters.
        """
        if isinstance(kmer, basestring):
            temp = kmer.encode('utf-8')
            return deref(self._ht_this).get_count(<char*>temp)
        # assume kmer is an integer representing the hash value
        else:
            return deref(self._ht_this).get_count(<uint64_t>kmer)

    def ksize(self):
        """k-mer size"""
        return deref(self._ht_this).ksize()

    def hashsizes(self):
        """Size of hash tables used."""
        return deref(self._ht_this).get_tablesizes()

    cdef _valid_sequence(self, sequence):
        """Validate sequence argument and convert it to bytes"""
        if len(sequence) < self.ksize():
            raise ValueError("sequence length ({}) must >= the hashtable "
                             "k-mer size ({})".format(len(sequence),
                                                      self.ksize()))
        return _bstring(sequence)

    def get_kmers(self, sequence):
        """Generate an ordered list of all k-mers in sequence."""
        data = self._valid_sequence(sequence)
        cdef vector[string] kmers
        deref(self._ht_this).get_kmers(data, kmers)
        return kmers

    def consume(self, sequence):
        """Increment the counts of all of the k-mers in the sequence."""
        data = self._valid_sequence(sequence)
        return deref(self._ht_this).consume_string(data)

    def get_kmer_counts(self, sequence):
        """Retrieve an ordered list of the counts of all k-mers in sequence."""
        data = self._valid_sequence(sequence)
        cdef vector[BoundedCounterType] counts
        deref(self._ht_this).get_kmer_counts(data, counts)
        return counts

    def get_min_count(self, sequence):
        """Get the smallest count of all the k-mers in the string."""
        data = self._valid_sequence(sequence)
        return deref(self._ht_this).get_min_count(data)

    def get_max_count(self, sequence):
        """Get the larget count of all the k-mers in the string."""
        data = self._valid_sequence(sequence)
        return deref(self._ht_this).get_max_count(data)

    def get_median_count(self, sequence):
        """median, average, and stddev of the k-mer counts in sequence."""
        data = self._valid_sequence(sequence)
        cdef BoundedCounterType med = 0
        cdef float average = 0
        cdef float stddev = 0

        deref(self._ht_this).get_median_count(data, med, average, stddev)
        return (med, average, stddev)

    def get_kmer_hashes(self, sequence):
        """Retrieve hashes of all k-mers in sequence.

        Hashes are returned in the same order as k-mers appear in sequence.
        """
        data = self._valid_sequence(sequence)
        cdef vector[HashIntoType] hashes
        deref(self._ht_this).get_kmer_hashes(data, hashes)
        return hashes

    def trim_on_abundance(self, sequence, abundance):
        """Trim sequence at first k-mer below the given abundance."""
        data = self._valid_sequence(sequence)
        trimmed_at = deref(self._ht_this).trim_on_abundance(data, abundance)
        return sequence[:trimmed_at], trimmed_at

    def trim_below_abundance(self, sequence, abundance):
        """Trim sequence at first k-mer above the given abundance."""
        data = self._valid_sequence(sequence)
        trimmed_at = deref(self._ht_this).trim_below_abundance(data, abundance)
        return sequence[:trimmed_at], trimmed_at

    def find_spectral_error_positions(self, sequence, max_count):
        """Identify positions of low-abundance k-mers."""
        data = self._valid_sequence(sequence)
        posns = (deref(self._ht_this).find_spectral_error_positions(data,
                                                                   max_count))
        return posns

    def consume_seqfile_with_reads_parser(self, read_parser):
        """Count all k-mers from read_parser."""
        cdef unsigned long long n_consumed = 0
        cdef unsigned int total_reads = 0

        cdef CPyReadParser_Object* parser = <CPyReadParser_Object*>read_parser

        deref(self._ht_this).consume_seqfile[CpFastxReader](parser.parser,
                                                           total_reads,
                                                           n_consumed)
        return total_reads, n_consumed

    def consume_seqfile(self, file_name):
        """Count all k-mers from file_name."""
        cdef unsigned long long n_consumed = 0
        cdef unsigned int total_reads = 0

        cdef FastxParserPtr parser = get_parser[CpFastxReader](_bstring(file_name))
        deref(self._ht_this).consume_seqfile[CpFastxReader](parser,
                                                           total_reads,
                                                           n_consumed)
        return total_reads, n_consumed

    def consume_seqfile_with_mask(self, file_name, Hashtable mask, int threshold=0):
        cdef unsigned long long n_consumed = 0
        cdef unsigned int total_reads = 0
        cdef FastxParserPtr parser = get_parser[CpFastxReader](_bstring(file_name))
        cdef CpHashtable * cmask = mask._ht_this.get()
        deref(self._ht_this).consume_seqfile_with_mask[CpFastxReader](parser,
                                                                     cmask,
                                                                     threshold,
                                                                     total_reads,
                                                                     n_consumed)
        return total_reads, n_consumed
                                                                     
    def consume_seqfile_banding(self, file_name, num_bands, band):
        """Count all k-mers from file_name."""
        cdef unsigned long long n_consumed = 0
        cdef unsigned int total_reads = 0
        cdef FastxParserPtr parser = get_parser[CpFastxReader](_bstring(file_name))
        deref(self._ht_this).consume_seqfile_banding[CpFastxReader](parser,
                                                                   num_bands,
                                                                   band,
                                                                   total_reads,
                                                                   n_consumed)
        return total_reads, n_consumed

    def consume_seqfile_banding_with_mask(self, file_name, num_bands, band,
                                          Hashtable mask, int threshold=0):
        cdef unsigned long long n_consumed = 0
        cdef unsigned int total_reads = 0
        cdef FastxParserPtr parser = get_parser[CpFastxReader](_bstring(file_name))
        cdef CpHashtable * cmask = mask._ht_this.get()
        deref(self._ht_this).consume_seqfile_banding_with_mask[CpFastxReader](parser,
                                                                     num_bands,
                                                                     band,
                                                                     cmask,
                                                                     threshold,
                                                                     total_reads,
                                                                     n_consumed)
        return total_reads, n_consumed

    def abundance_distribution(self, file_name, Hashtable tracking):
        """Calculate the k-mer abundance distribution over reads in file_name."""
        cdef FastxParserPtr parser = get_parser[CpFastxReader](_bstring(file_name))
        cdef CpHashtable * cptracking = tracking._ht_this.get()
        cdef uint64_t * x = deref(self._ht_this).\
                abundance_distribution[CpFastxReader](parser, cptracking)
        abunds = []
        for i in range(MAX_BIGCOUNT):
            abunds.append(x[i])
        return abunds

    def abundance_distribution_with_reads_parser(self, object read_parser, Hashtable tracking):
        """Calculate the k-mer abundance distribution over reads."""

        cdef CpHashtable * cptracking = tracking._ht_this.get()
 
        cdef CPyReadParser_Object* parser
        parser = <CPyReadParser_Object*>read_parser
        cdef uint64_t * x = deref(self._ht_this).abundance_distribution[CpFastxReader](
                parser.parser, cptracking)
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
    def __cinit__(self, int k, int starting_size):
        # starting size has to be a power of two
        power_of_two = ((starting_size & (starting_size - 1) == 0) and
                        (starting_size != 0))
        if not power_of_two:
            raise ValueError("starting_size has to be a power of two.")
        if type(self) is QFCounttable:
            self._qf_this = make_shared[CpQFCounttable](k, <uint64_t>log(starting_size, 2))
            self._ht_this = <shared_ptr[CpHashtable]>self._qf_this

    @classmethod
    def load(cls, file_name):
        """Load the graph from the specified file."""
        cdef QFCounttable table = cls(1, 1)
        deref(table._qf_this).load(_bstring(file_name))
        return table

cdef class Counttable(Hashtable):

    def __cinit__(self, int k, int starting_size, int n_tables):
        cdef vector[uint64_t] primes
        if type(self) is Counttable:
            primes = get_n_primes_near_x(n_tables, starting_size)
            self._ct_this = make_shared[CpCounttable](k, primes)
            self._ht_this = <shared_ptr[CpHashtable]>self._ct_this


cdef class SmallCounttable(Hashtable):

    def __cinit__(self, int k, int starting_size, int n_tables):
        cdef vector[uint64_t] primes
        if type(self) is SmallCounttable:
            primes = get_n_primes_near_x(n_tables, starting_size)
            self._st_this = make_shared[CpSmallCounttable](k, primes)
            self._ht_this = <shared_ptr[CpHashtable]>self._st_this

    def get_raw_tables(self):
        cdef uint8_t ** table_ptrs = deref(self._st_this).get_raw_tables()
        cdef vector[uint64_t] sizes = deref(self._st_this).get_tablesizes()
        for i in range(len(sizes)):
            sizes[i] = sizes[i] / 2 + 1
        return self._get_raw_tables(table_ptrs, sizes)


cdef class Nodetable(Hashtable):

    def __cinit__(self, int k, int starting_size, int n_tables):
        cdef vector[uint64_t] primes
        if type(self) is Nodetable:
            primes = get_n_primes_near_x(n_tables, starting_size)
            self._nt_this = make_shared[CpNodetable](k, primes)
            self._ht_this = <shared_ptr[CpHashtable]>self._nt_this


cdef class Hashgraph(Hashtable):

    def neighbors(self, str kmer):
        '''Get a list of neighbor nodes for this k-mer.'''
        cdef Traverser traverser = Traverser(self)
        return list(traverser.neighbors(kmer))

    def calc_connected_graph_size(self, str kmer, max_size=0,
                                  break_on_circumference=False):
        '''Find the number of nodes connected to this k-mer.'''
        cdef CpKmer start = deref(self._hg_this).build_kmer(_bstring(kmer))
        cdef uint64_t _size = 0
        cdef uint64_t _max_size = max_size
        cdef bool _break = break_on_circumference
        cdef KmerSet keeper
        cdef CpHashgraph * ptr = self._hg_this.get() # need tmp ref for nogil

        with nogil:
            deref(ptr).calc_connected_graph_size(start, _size,
                                                 keeper, _max_size,
                                                 _break)
        return _size

    def kmer_degree(self, str kmer):
        '''Calculate the number of immediate neighbors this k-mer has
        the graph.'''
        return deref(self._hg_this).kmer_degree(_bstring(kmer))

    def count_kmers_within_radius(self, str kmer, int radius, int max_count=0):
        '''Calculate the number of neighbors with given radius in the graph.'''
        cdef unsigned int n
        cdef uint32_t _radius = radius
        cdef uint32_t _max_count = max_count
        cdef CpKmer start = deref(self._hg_this).build_kmer(_bstring(kmer))
        cdef KmerSet seen
        cdef CpHashgraph * ptr = self._hg_this.get()
        with nogil:
            n = deref(ptr).traverse_from_kmer(start, _radius,
                                              seen, _max_count)
        return n

    def find_high_degree_nodes(self, str sequence):
        '''Examine the given sequence for degree > 2 nodes and add to
        list; used in graph contraction.'''
        cdef HashSet hdns = HashSet(self.ksize())
        data = self._valid_sequence(sequence)
        deref(self._hg_this).find_high_degree_nodes(data, 
                                                   hdns.hs)
        return hdns


    def traverse_linear_path(self, str kmer, HashSet hdns, 
                             Nodegraph stop_filter=None):
        '''Traverse the path through the graph starting with the given
        k-mer and avoiding high-degree nodes, finding (and returning)
        traversed k-mers and any encountered high-degree nodes.'''
        cdef set[HashIntoType] adj
        cdef set[HashIntoType] visited
        cdef CpKmer cpkmer = CpKmer(_bstring(kmer), self.ksize())
        cdef CpNodegraph * _stop_filter = stop_filter._ng_this.get()
        cdef int size = deref(self._hg_this).traverse_linear_path(cpkmer,
                                                                 adj,
                                                                 visited,
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
            
    def find_all_tags_list(self, str kmer):
        '''Find all tags within range of the given k-mer, return as list'''
        if len(kmer) != self.ksize():
            raise ValueError("k-mer length must equal the counting "\
                             "table k-mer size")
        cdef HashSet result = HashSet(self.ksize())
        cdef CpKmer start = deref(self._hg_this).build_kmer(_bstring(kmer))

        with nogil:
            # partition->find_all_tags(start_kmer, result.hs, all_tags)
            pass

        return result


    def consume_seqfile_and_tag(self, str filename):
        '''Consume all sequences in a FASTA/FASTQ file and tag the resulting
        graph.'''
        cdef unsigned long long n_consumed = 0
        cdef unsigned int total_reads = 0
        cdef string _filename = _bstring(filename)

        deref(self._hg_this).consume_seqfile_and_tag[CpFastxReader](_filename,
                                                                   total_reads,
                                                                   n_consumed)
        return total_reads, n_consumed
    
    def print_tagset(self, str filename):
        '''Print out all of the tags.'''
        deref(self._hg_this).print_tagset(_bstring(filename))
    
    def add_tag(self, object kmer):
        '''Add a k-mer to the tagset.'''
        if isinstance(kmer, basestring):
            deref(self._hg_this).add_tag(deref(self._hg_this).hash_dna(_bstring(kmer)))
        else:
            return deref(self._hg_this).add_tag(<uint64_t>kmer)
    
    def get_tagset(self):
        '''Get all tagged k-mers as DNA strings.'''
        cdef HashIntoType st
        cdef list all_tags = []
        for st in deref(self._hg_this).all_tags:
            all_tags.append(deref(self._hg_this).unhash_dna(st))
        return all_tags

    def iter_tagset(self):
        '''Get all tagged k-mers as DNA strings.'''
        cdef HashIntoType st
        for st in deref(self._hg_this).all_tags:
            yield deref(self._hg_this).unhash_dna(st)

    def load_tagset(self, str filename, clear_tags=False):
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

    def do_subset_partition(self):
        '''Partition the graph starting from a given subset of tags.'''
        pass
    
    def find_all_tags(self):
        '''Starting from the given k-mer, find all closely connected tags.'''
        pass
    
    def assign_partition_id(self):
        '''Assign a partition ID to a given tag.'''
        pass
    
    def output_partitions(self):
        '''Write out sequences in given filename to another file, annotating '''
        '''with partition IDs.'''
        pass
    
    def load_partitionmap(self):
        '''Load a partitionmap for a given subset.'''
        pass

    def save_partitionmap(self):
        '''Save a partitionmap for the given subset.'''
        pass
    
    def _validate_partitionmap(self):
        '''Run internal validation checks.'''
        pass
    
    def consume_seqfile_and_tag_with_reads_parser(self, object read_parser):
        '''Count all k-mers using the given reads parser'''
        cdef unsigned long long n_consumed = 0
        cdef unsigned int total_reads = 0
        cdef CPyReadParser_Object * parser_o = <CPyReadParser_Object*>read_parser
        cdef FastxParserPtr parser = parser_o.parser
        cdef CpHashgraph * ptr = self._hg_this.get()

        deref(ptr).consume_seqfile_and_tag_readparser[CpFastxReader](parser,
                                                            total_reads,
                                                            n_consumed)
        return total_reads, n_consumed
    
    def consume_partitioned_fasta(self, filename):
        '''Count all k-mers in a given file'''
        cdef unsigned long long n_consumed = 0
        cdef unsigned int total_reads = 0
        cdef string _filename = _bstring(filename)
        deref(self._hg_this).consume_partitioned_fasta[CpFastxReader](_filename,
                                                                     total_reads,
                                                                     n_consumed)
        return total_reads, n_consumed
    
    def merge_subset(self):
        '''Merge the given subset into this one.'''
        pass
    
    def merge_subset_from_disk(self):
        '''Merge the given subset (filename) into this one.'''
        pass
    
    def count_partitions(self):
        '''Count the number of partitions in the master partitionmap.'''
        pass
    
    def subset_count_partitions(self):
        '''Count the number of partitions in this subset partitionmap.'''
        pass

    def subset_partition_size_distribution(self):
        '''Get the size distribution of partitions in this subset.'''
        pass

    def save_subset_partitionmap(self):
        '''Save the partition map for this subset.'''
        pass

    def load_subset_partitionmap(self):
        '''Save the partition map for this subset.'''
        pass
    
    def _validate_subset_partitionmap(self):
        '''Run internal validation checks on this subset.'''
        pass
    
    def set_partition_id(self):
        '''Set the partition ID for this tag.'''
        pass

    def join_partitions(self):
        '''Join the partitions of these two tags.'''
        pass
    
    def get_partition_id(self):
        '''Get the partition ID of this tag.'''
        pass
    
    def repartition_largest_partition(self):
        '''Repartition the largest partition (in the face of stop tags).'''
        pass

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
        return sequence[:trim_at]

    def add_stop_tag(self, object kmer):
        '''Add this k-mer as a stop tag.'''
        if isinstance(kmer, basestring):
            deref(self._hg_this).add_stop_tag(deref(self._hg_this).hash_dna(_bstring(kmer)))
        else:
            return deref(self._hg_this).add_stop_tag(<uint64_t>kmer)
    
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

    def __cinit__(self, int k, int starting_size, int n_tables,
                  primes=[]):
        cdef vector[uint64_t] _primes
        if type(self) is Countgraph:
            if primes:
                _primes = primes
            else:
                _primes = get_n_primes_near_x(n_tables, starting_size)
            self._cg_this = make_shared[CpCountgraph](k, _primes)
            self._hg_this = <shared_ptr[CpHashgraph]>self._cg_this
            self._ht_this = <shared_ptr[CpHashtable]>self._hg_this


cdef class SmallCountgraph(Hashgraph):

    def __cinit__(self, int k, int starting_size, int n_tables,
                  primes=[]):
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

    def __cinit__(self, int k, int starting_size, int n_tables,
                  primes=[]):
        cdef vector[uint64_t] _primes
        if type(self) is Nodegraph:
            if primes:
                _primes = primes
            else:
                _primes = get_n_primes_near_x(n_tables, starting_size)
            self._ng_this = make_shared[CpNodegraph](k, _primes)
            self._hg_this = <shared_ptr[CpHashgraph]>self._ng_this
            self._ht_this = <shared_ptr[CpHashtable]>self._hg_this

    def update_from(self, Nodegraph other):
        deref(self._ng_this).update_from(deref(other._ng_this))
