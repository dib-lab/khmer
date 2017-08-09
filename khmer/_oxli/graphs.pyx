# cython: c_string_type=unicode, c_string_encoding=utf8
from math import log

from cython.operator cimport dereference as deref
from libc.stdint cimport uint64_t

from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector

from utils cimport _bstring
from utils import get_n_primes_near_x
from parsing cimport (CpFastxReader, CPyReadParser_Object, get_parser,
                      CpReadParser, FastxParserPtr)
from oxli_types cimport MAX_BIGCOUNT
from .._khmer import Countgraph as PyCountgraph
from .._khmer import Nodegraph as PyNodegraph
from .._khmer import GraphLabels as PyGraphLabels
from .._khmer import ReadParser


CYTHON_TABLES = (Hashtable, Nodetable, Counttable, SmallCounttable,
                 QFCounttable)
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


cdef CpHashtable * hashtable_arg_shim(object table,
                                      allowed=(PyNodegraph, PyCountgraph,
                                               Nodetable, Counttable,
                                               SmallCounttable, QFCounttable)):
    cdef CPyHashtable_Object* cpyhashtable
    cdef CpHashtable * hashtable
                                          
    if isinstance(table, allowed):
        if isinstance(table, CYTHON_TABLES):
            hashtable = (<Hashtable>table).c_table.get()
        else:
            cpyhashtable = <CPyHashtable_Object*>table
            hashtable = cpyhashtable.hashtable
    else:
        raise ValueError('Expected one of {0}, '\
                         'got {1} instead.'.format(allowed, type(table)))

    return hashtable


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
            return deref(self.c_table).add(<char*>temp)
        # assume kmer is an integer representing the hash value
        else:
            return deref(self.c_table).add(<uint64_t>kmer)

    def hash(self, kmer):
        """Compute the hash of this k-mer."""
        if len(kmer) != self.ksize():
            raise ValueError("Expected k-mer length {}"
                             " but got {}.".format(self.ksize(), len(kmer)))
        data = _bstring(kmer)
        return deref(self.c_table).hash_dna(data)

    def reverse_hash(self, kmer_hash):
        """Turn a k-mer hash back into a DNA k-mer, if possible."""
        return deref(self.c_table).unhash_dna(kmer_hash)

    def get(self, kmer):
        """Retrieve the count for the given k-mer.

        `kmer` can be either a string or an integer representing the hashed
        value of the kmer.

        For Nodetables and Counttables, this function will fail if the
        supplied k-mer contains non-ACGT characters.
        """
        if isinstance(kmer, basestring):
            temp = kmer.encode('utf-8')
            return deref(self.c_table).get_count(<char*>temp)
        # assume kmer is an integer representing the hash value
        else:
            return deref(self.c_table).get_count(<uint64_t>kmer)

    def ksize(self):
        """k-mer size"""
        return deref(self.c_table).ksize()

    def hashsizes(self):
        """Size of hash tables used."""
        return deref(self.c_table).get_tablesizes()

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
        deref(self.c_table).get_kmers(data, kmers)
        return kmers

    def consume(self, sequence):
        """Increment the counts of all of the k-mers in the sequence."""
        data = self._valid_sequence(sequence)
        return deref(self.c_table).consume_string(data)

    def get_kmer_counts(self, sequence):
        """Retrieve an ordered list of the counts of all k-mers in sequence."""
        data = self._valid_sequence(sequence)
        cdef vector[BoundedCounterType] counts
        deref(self.c_table).get_kmer_counts(data, counts)
        return counts

    def get_min_count(self, sequence):
        """Get the smallest count of all the k-mers in the string."""
        data = self._valid_sequence(sequence)
        return deref(self.c_table).get_min_count(data)

    def get_max_count(self, sequence):
        """Get the larget count of all the k-mers in the string."""
        data = self._valid_sequence(sequence)
        return deref(self.c_table).get_max_count(data)

    def get_median_count(self, sequence):
        """median, average, and stddev of the k-mer counts in sequence."""
        data = self._valid_sequence(sequence)
        cdef BoundedCounterType med = 0
        cdef float average = 0
        cdef float stddev = 0

        deref(self.c_table).get_median_count(data, med, average, stddev)
        return (med, average, stddev)

    def get_kmer_hashes(self, sequence):
        """Retrieve hashes of all k-mers in sequence.

        Hashes are returned in the same order as k-mers appear in sequence.
        """
        data = self._valid_sequence(sequence)
        cdef vector[HashIntoType] hashes
        deref(self.c_table).get_kmer_hashes(data, hashes)
        return hashes

    def trim_on_abundance(self, sequence, abundance):
        """Trim sequence at first k-mer below the given abundance."""
        data = self._valid_sequence(sequence)
        trimmed_at = deref(self.c_table).trim_on_abundance(data, abundance)
        return sequence[:trimmed_at], trimmed_at

    def trim_below_abundance(self, sequence, abundance):
        """Trim sequence at first k-mer above the given abundance."""
        data = self._valid_sequence(sequence)
        trimmed_at = deref(self.c_table).trim_below_abundance(data, abundance)
        return sequence[:trimmed_at], trimmed_at

    def find_spectral_error_positions(self, sequence, max_count):
        """Identify positions of low-abundance k-mers."""
        data = self._valid_sequence(sequence)
        posns = (deref(self.c_table).find_spectral_error_positions(data,
                                                                   max_count))
        return posns

    def consume_seqfile_with_reads_parser(self, read_parser):
        """Count all k-mers from read_parser."""
        cdef unsigned long long n_consumed = 0
        cdef unsigned int total_reads = 0

        cdef CPyReadParser_Object* parser = <CPyReadParser_Object*>read_parser

        deref(self.c_table).consume_seqfile[CpFastxReader](parser.parser,
                                                           total_reads,
                                                           n_consumed)
        return total_reads, n_consumed

    def consume_seqfile(self, file_name):
        """Count all k-mers from file_name."""
        cdef unsigned long long n_consumed = 0
        cdef unsigned int total_reads = 0

        cdef FastxParserPtr parser = get_parser[CpFastxReader](_bstring(file_name))
        deref(self.c_table).consume_seqfile[CpFastxReader](parser,
                                                           total_reads,
                                                           n_consumed)
        return total_reads, n_consumed

    def consume_seqfile_with_mask(self, file_name, mask, threshold=0):
        cdef unsigned long long n_consumed = 0
        cdef unsigned int total_reads = 0
        cdef FastxParserPtr parser = get_parser[CpFastxReader](_bstring(file_name))
        cdef CpHashtable * cmask = hashtable_arg_shim(mask)
        deref(self.c_table).consume_seqfile_with_mask[CpFastxReader](parser,
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
        deref(self.c_table).consume_seqfile_banding[CpFastxReader](parser,
                                                                   num_bands,
                                                                   band,
                                                                   total_reads,
                                                                   n_consumed)
        return total_reads, n_consumed

    def consume_seqfile_banding_with_mask(self, file_name, num_bands, band,
                                          mask, threshold=0):
        cdef unsigned long long n_consumed = 0
        cdef unsigned int total_reads = 0
        cdef FastxParserPtr parser = get_parser[CpFastxReader](_bstring(file_name))
        cdef CpHashtable * cmask = hashtable_arg_shim(mask)
        deref(self.c_table).consume_seqfile_banding_with_mask[CpFastxReader](parser,
                                                                     num_bands,
                                                                     band,
                                                                     cmask,
                                                                     threshold,
                                                                     total_reads,
                                                                     n_consumed)
        return total_reads, n_consumed

    def abundance_distribution(self, file_name, tracking):
        """Calculate the k-mer abundance distribution over reads in file_name."""
        cdef FastxParserPtr parser = get_parser[CpFastxReader](_bstring(file_name))
        cdef CpHashtable * cptracking = hashtable_arg_shim(tracking,
                                                      allowed=(PyNodegraph, Nodetable))
        cdef uint64_t * x = deref(self.c_table).\
                abundance_distribution[CpFastxReader](parser, cptracking)
        abunds = []
        for i in range(MAX_BIGCOUNT):
            abunds.append(x[i])
        return abunds

    def abundance_distribution_with_reads_parser(self, read_parser, tracking):
        """Calculate the k-mer abundance distribution over reads."""
        cdef CpHashtable * cptracking = hashtable_arg_shim(tracking,
                                                      allowed=(PyNodegraph, Nodetable))
 
        cdef CPyReadParser_Object* parser
        parser = <CPyReadParser_Object*>read_parser
        cdef uint64_t * x = deref(self.c_table).abundance_distribution[CpFastxReader](
                parser.parser, cptracking)
        abunds = []
        for i in range(MAX_BIGCOUNT):
            abunds.append(x[i])
        return abunds

    def save(self, file_name):
        """Save the graph to the specified file."""
        deref(self.c_table).save(_bstring(file_name))

    @classmethod
    def load(cls, file_name):
        """Load the graph from the specified file."""
        cdef Hashtable table = cls(1, 1, 1)
        deref(table.c_table).load(_bstring(file_name))
        return table

    def n_unique_kmers(self):
        """Estimate of the number of unique kmers stored."""
        return deref(self.c_table).n_unique_kmers()

    def n_occupied(self):
        """Estimate of the number of occupied slots in the storage."""
        return deref(self.c_table).n_occupied()

    def n_tables(self):
        """Number of tables used in the storage."""
        return deref(self.c_table).n_tables()

    def set_use_bigcount(self, bigcount):
        deref(self.c_table).set_use_bigcount(bigcount)

    def get_use_bigcount(self):
        return deref(self.c_table).get_use_bigcount()


cdef class QFCounttable(Hashtable):
    def __cinit__(self, int k, int starting_size):
        # starting size has to be a power of two
        power_of_two = ((starting_size & (starting_size - 1) == 0) and
                        (starting_size != 0))
        if not power_of_two:
            raise ValueError("starting_size has to be a power of two.")
        if type(self) is QFCounttable:
            self.c_table.reset(<CpHashtable*>new CpQFCounttable(k, int(log(starting_size, 2))))

    @classmethod
    def load(cls, file_name):
        """Load the graph from the specified file."""
        cdef Hashtable table = cls(1, 1)
        deref(table.c_table).load(_bstring(file_name))
        return table

cdef class Counttable(Hashtable):

    def __cinit__(self, int k, int starting_size, int n_tables):
        if type(self) is Counttable:
            primes = get_n_primes_near_x(n_tables, starting_size)
            self.c_table.reset(<CpHashtable*>new CpCounttable(k, primes))


cdef class SmallCounttable(Hashtable):

    def __cinit__(self, int k, int starting_size, int n_tables):
        if type(self) is SmallCounttable:
            primes = get_n_primes_near_x(n_tables, starting_size)
            self.c_table.reset(<CpHashtable*>new CpSmallCounttable(k, primes))


cdef class Nodetable(Hashtable):

    def __cinit__(self, int k, int starting_size, int n_tables):
        if type(self) is Nodetable:
            primes = get_n_primes_near_x(n_tables, starting_size)
            self.c_table.reset(<CpHashtable*>new CpNodetable(k, primes))


cdef class Hashgraph(Hashtable):

    def neighbors(self):
        '''Get a list of neighbor nodes for this k-mer.'''
        pass

    def calc_connected_graph_size(self):
        '''Find the number of nodes connected to this k-mer.'''
        pass

    def kmer_degree(self):
        '''Calculate the number of immediate neighbors this k-mer has
        the graph.'''
        pass

    def count_kmers_within_radius(self):
        '''Calculate the number of neighbors with given radius in the graph.'''
        pass

    def find_high_degree_nodes(self, kmer):
        '''Examine the given sequence for degree > 2 nodes and add to
        list; used in graph contraction.'''
        pass

    def traverse_linear_path(self):
        '''Traverse the path through the graph starting with the given
        "k-mer and avoiding high-degree nodes, finding (and returning)
        "traversed k-mers and any encountered high-degree nodes.'''
        pass

    def assemble_linear_path(self):
        '''Assemble a purely linear path starting with the given
        "k-mer, returning traversed k-mers and any encountered high-degree
        "nodes.'''
        pass


    def consume_and_tag(self):
        "Consume a sequence and tag it."
        pass

    def get_tags_and_positions(self):
        "Retrieve tags and their positions in a sequence."
        pass
    
    def find_all_tags_list(self):
        "Find all tags within range of the given k-mer, return as list"
        pass

    def consume_seqfile_and_tag(self):
        "Consume all sequences in a FASTA/FASTQ file and tag the resulting "
        "graph."
        pass
    
    def print_tagset(self):
        "Print out all of the tags."
        pass
    
    def add_tag(self):
        "Add a k-mer to the tagset."
        pass
    
    def get_tagset(self):
        "Get all tagged k-mers as DNA strings."
        pass

    def load_tagset(self):
        "Load tags from a file."
        pass
    
    def save_tagset(self):
        "Save tags to a file."
        pass
    
    def n_tags(self):
        "Return the count of all tags."
        pass
    
    def divide_tags_into_subsets(self):
        "Divide tags equally up into subsets of given size."
        pass
    
    @property
    def tag_density(self):
        "Get the tagging density."
        pass
    
    def _set_tag_density(self):
        "Set the tagging density."
        pass

    def do_subset_partition(self):
        "Partition the graph starting from a given subset of tags."
        pass
    
    def find_all_tags(self):
        "Starting from the given k-mer, find all closely connected tags."
        pass
    
    def assign_partition_id(self):
        "Assign a partition ID to a given tag."
        pass
    
    def output_partitions(self):
        "Write out sequences in given filename to another file, annotating "
        "with partition IDs."
        pass
    
    def load_partitionmap(self):
        "Load a partitionmap for a given subset."
        pass

    def save_partitionmap(self):
        "Save a partitionmap for the given subset."
        pass
    
    def _validate_partitionmap(self):
        "Run internal validation checks."
        pass
    
    def consume_seqfile_and_tag_with_reads_parser(self):
        "Count all k-mers using the given reads parser"
        pass
    
    def consume_partitioned_fasta(self):
        "Count all k-mers in a given file"
        pass
    
    def merge_subset(self):
        "Merge the given subset into this one."
        pass
    
    def merge_subset_from_disk(self):
        "Merge the given subset (filename) into this one."
        pass
    
    def count_partitions(self):
        "Count the number of partitions in the master partitionmap."
        pass
    
    def subset_count_partitions(self):
        "Count the number of partitions in this subset partitionmap."
        pass

    def subset_partition_size_distribution(self):
        "Get the size distribution of partitions in this subset."
        pass

    def save_subset_partitionmap(self):
        "Save the partition map for this subset."
        pass

    def load_subset_partitionmap(self):
        "Save the partition map for this subset."
        pass
    
    def _validate_subset_partitionmap(self):
        "Run internal validation checks on this subset."
        pass
    
    def set_partition_id(self):
        "Set the partition ID for this tag."
        pass

    def join_partitions(self):
        "Join the partitions of these two tags."
        pass
    
    def get_partition_id(self):
        "Get the partition ID of this tag."
        pass
    
    def repartition_largest_partition(self):
        "Repartition the largest partition (in the face of stop tags)."
        pass

    def load_stop_tags(self):
        "Load the set of stop tags."
        
    def save_stop_tags(self):
        "Save the set of stop tags."
        pass

    def print_stop_tags(self):
        "Print out the set of stop tags."
        pass
    
    def trim_on_stoptags(self):
        "Trim the reads on the given stop tags."
        pass

    def add_stop_tag(self):
        "Add this k-mer as a stop tag."
        pass

    def get_stop_tags(self):
        "Return a DNA list of all of the stop tags."
        pass
