# cython: c_string_type=unicode, c_string_encoding=utf8
from math import log2

from cython.operator cimport dereference as deref
from libc.stdint cimport uint64_t

from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector

from .utils cimport _bstring
from graphs cimport CpQFCounttable
from parsing cimport CpFastxReader, CPyReadParser_Object
from oxli_types cimport MAX_BIGCOUNT
from .._khmer import Countgraph, Nodegraph, GraphLabels, Nodetable, ReadParser


cdef CpHashgraph * get_hashgraph_ptr(object graph):
    if not (isinstance(graph, Countgraph) or isinstance(graph, Nodegraph)):
        return NULL

    cdef CPyHashgraph_Object* ptr = <CPyHashgraph_Object*> graph
    return deref(ptr).hashgraph


cdef CpLabelHash * get_labelhash_ptr(object labels):
    if not isinstance(labels, GraphLabels):
        return NULL

    cdef CPyGraphLabels_Object * ptr = <CPyGraphLabels_Object*> labels
    return deref(ptr).labelhash


cdef class QFCounttable:
    cdef unique_ptr[CpQFCounttable] c_table
    def __cinit__(self, int k, int starting_size):
        # starting size has to be a power of two
        power_of_two = ((starting_size & (starting_size - 1) == 0) and
                        (starting_size != 0))
        if not power_of_two:
            raise ValueError("starting_size has to be a power of two.")
        self.c_table.reset(new CpQFCounttable(k, int(log2(starting_size))))

    def add(self, kmer):
        """Increment the count of this k-mer.

        Synonym for 'count'.
        """
        return self.count(kmer)

    def count(self, kmer):
        """Increment the count of this k-mer

        `kmer` can be either a string or an integer representing the hashed
        value of the kmer.
        """
        if isinstance(kmer, str):
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
        if isinstance(kmer, str):
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

    def consume_seqfile(self, file_name):
        """Count all k-mers from file_name."""
        cdef unsigned long long n_consumed = 0
        cdef unsigned int total_reads = 0

        read_parser = ReadParser(file_name)
        cdef CPyReadParser_Object* parser = <CPyReadParser_Object*>read_parser
        deref(self.c_table).consume_seqfile[CpFastxReader](parser.parser,
                                                           total_reads,
                                                           n_consumed)

    def abundance_distribution(self, file_name, tracking):
        """Calculate the k-mer abundance distribution over reads in file_name."""
        cdef CPyReadParser_Object* parser
        if isinstance(file_name, str):
            read_parser = ReadParser(file_name)
            parser = <CPyReadParser_Object*>read_parser

        else:
            raise ValueError('Expected file_name to be string, '
                             'got {} instead.'.format(type(file_name)))

        cdef CPyHashtable_Object* hashtable
        if isinstance(tracking, (Nodetable, Nodegraph)):
            hashtable = <CPyHashtable_Object*>tracking
        else:
            raise ValueError('Expected `tracking` to be a Nodetable or '
                             'Nodegraph, got {} instead.'.format(type(tracking)))

        cdef uint64_t * x = deref(self.c_table).abundance_distribution[CpFastxReader](
                parser.parser, hashtable.hashtable)
        abunds = []
        for i in range(MAX_BIGCOUNT):
            abunds.append(x[i])
        return abunds

    def abundance_distribution_with_reads_parser(self, read_parser, tracking):
        """Calculate the k-mer abundance distribution over reads."""
        cdef CPyHashtable_Object* hashtable
        if isinstance(tracking, (Nodetable, Nodegraph)):
            hashtable = <CPyHashtable_Object*>tracking
        else:
            raise ValueError('Expected `tracking` to be a Nodetable or '
                             'Nodegraph, got {} instead.'.format(type(tracking)))

        cdef CPyReadParser_Object* parser
        parser = <CPyReadParser_Object*>read_parser
        cdef uint64_t * x = deref(self.c_table).abundance_distribution[CpFastxReader](
                parser.parser, hashtable.hashtable)
        abunds = []
        for i in range(MAX_BIGCOUNT):
            abunds.append(x[i])
        return abunds
