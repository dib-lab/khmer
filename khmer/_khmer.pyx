'''
This file is part of khmer, http://github.com/ged-lab/khmer/, and is
Copyright (C) Michigan State University, 2009-2014. It is licensed under
the three-clause BSD license; see doc/LICENSE.txt.
Contact: khmer-project@idyll.org
'''

from cython.operator cimport dereference as deref
from libc.limits cimport UINT_MAX

# FIXME: ugly, but Cython has no way to check #define macros
cdef int MAX_BIGCOUNT_C = 65535
cdef int MAX_COUNT_C = 255


cdef class _LabelHash:
    pass


cdef class Hashbits:
    cdef CppHashbits *thisptr

    def __cinit__(self, WordLength k, vector[unsigned long long int]& sizes):
        self.thisptr = new CppHashbits(k, sizes)

    def __dealloc__(self):
        del self.thisptr

    def save(self, string filename):
        self.thisptr.save(filename)

    def consume_fasta(self, const string filename, callback_obj=None):
        cdef unsigned long long n_consumed = 0
        cdef unsigned int total_reads = 0
        cdef CallbackFn _report_fn = NULL

#        try {
        self.thisptr.consume_fasta(filename, total_reads, n_consumed,
                                   _report_fn, <void*>callback_obj)
#        } catch (_khmer_signal &e) {
#            return NULL;
#        }

        return total_reads, n_consumed

    def calc_connected_graph_size(self, const char *_kmer, unsigned int max_size=0,
                                        bool break_on_circum=False):
        cdef unsigned long long size = 0
        cdef SeenSet keeper
        self.thisptr.calc_connected_graph_size(_kmer, size, keeper, max_size,
                                                break_on_circum)
        return size

    def consume(self, const string& long_str):
        return self.thisptr.consume_string(long_str)

    def consume_fasta_and_tag(self, const string filename, callback_obj=None):
        cdef unsigned long long n_consumed = 0
        cdef unsigned int total_reads = 0
        cdef CallbackFn _report_fn = NULL

#        try {
        self.thisptr.consume_fasta_and_tag(filename, total_reads,
                                           n_consumed, _report_fn,
                                           <void*>callback_obj)
#        } catch (_khmer_signal &e) {
#            return NULL;
#        }

        return total_reads, n_consumed

    def do_subset_partition(self, HashIntoType start_kmer=0,
                            HashIntoType end_kmer=0,
                            bool break_on_stop_tags=False,
                            bool stop_big_traversals=False,
                            callback_obj=None):
        cdef CallbackFn _report_fn = NULL
        #try {
        cdef SubsetPartition *subset_p = new SubsetPartition(self.thisptr)
        subset_p.do_partition(start_kmer, end_kmer, break_on_stop_tags,
                               stop_big_traversals,
                               _report_fn, <void*>callback_obj)
        #} catch (_khmer_signal &e) {
        #    return NULL;
        #}

        return spi_factory(subset_p)

    def subset_count_partitions(self, SubsetPartitionInfo subset_p):
        cdef unsigned int n_partitions = 0, n_unassigned = 0
        subset_p.thisptr.count_partitions(n_partitions, n_unassigned)

        return n_partitions, n_unassigned

    def output_partitions(self, const string filename,
                          const string output_file,
                          bool output_unassigned=False,
                          callback_obj=None):
        cdef unsigned int n_partitions = 0
        cdef CallbackFn _report_fn = NULL

        #try {
        cdef SubsetPartition * subset_p = self.thisptr.partition
        n_partitions = subset_p.output_partitioned_file(filename,
                            output_file,
                            output_unassigned,
                            _report_fn,
                            <void*>callback_obj)
        #} catch (_khmer_signal &e) {
        #    return NULL;
        #}

        return n_partitions

    def merge_subset(self, SubsetPartitionInfo subset_p):
        self.thisptr.partition.merge(subset_p.thisptr)

    def find_all_tags(self, const char *kmer_s):
        if len(kmer_s) < self.thisptr.ksize():
            raise ValueError("starting kmer is smaller than the K size of the hashbits")

        cdef pre_partition_info * ppi

        cdef HashIntoType kmer, kmer_f, kmer_r
        kmer = _hash(kmer_s, self.thisptr.ksize(), kmer_f, kmer_r)

        ppi = new pre_partition_info(kmer)
        self.thisptr.partition.find_all_tags(kmer_f, kmer_r, ppi.tagged_kmers,
                                             self.thisptr.all_tags, False, False)
        self.thisptr.add_kmer_to_tags(kmer)

        return ppi_factory(ppi)

    def assign_partition_id(self, PrePartitionInfo ppi):
        return self.thisptr.partition.assign_partition_id(
                    ppi.thisptr.kmer, ppi.thisptr.tagged_kmers)

    def _get_tag_density(self):
        return self.thisptr._get_tag_density()

    def _set_tag_density(self, unsigned int d):
        self.thisptr._set_tag_density(d)

    def n_occupied(self, HashIntoType start=0, HashIntoType stop=0):
        return self.thisptr.n_occupied(start, stop)

    def get(self, val):
        if isinstance(val, int) or isinstance(val, long):
            return self.thisptr.get_count(<HashIntoType>val);
        else:
            return self.thisptr.get_count(<char *>val);

    def n_unique_kmers(self,  HashIntoType start=0, HashIntoType stop=0):
        return self.thisptr.n_kmers(start, stop)

    def count(self, const char* kmer):
        """Count the given kmer"""
        if len(kmer) != self.thisptr.ksize():
            raise ValueError("k-mer length must be the same as the hashtable k-size")
        self.thisptr.count(kmer)
        return 1

    def filter_if_present(self, const string filename, const string output,
                                callback_obj=None):
        cdef CallbackFn _report_fn = NULL
        #try {
        self.thisptr.filter_if_present(filename, output, _report_fn, <void*>callback_obj)
        #} catch (_khmer_signal &e) {
        #    return NULL;
        #}


    def consume_partitioned_fasta(self, const string filename, callback_obj=None):
        cdef unsigned long long n_consumed
        cdef unsigned int total_reads
        cdef CallbackFn _report_fn = NULL

        #try {
        self.thisptr.consume_partitioned_fasta(filename, total_reads, n_consumed,
                                                _report_fn, <void*>callback_obj)
        #} catch (_khmer_signal) {
        #    return NULL;
        #}

        return total_reads, n_consumed

    def count_partitions(self):
        cdef unsigned int n_partitions = 0, n_unassigned = 0
        self.thisptr.partition.count_partitions(n_partitions, n_unassigned)

        return n_partitions, n_unassigned

    def get_partition_id(self, string kmer):
        return self.thisptr.partition.get_partition_id(kmer)

    def join_partitions(self, PartitionID p1, PartitionID p2):
        return self.thisptr.partition.join_partitions(p1, p2)


    def count_kmers_within_radius(self, const char *kmer, unsigned int radius,
                                        unsigned int max_count=0):
        cdef HashIntoType kmer_f, kmer_r

        _hash(kmer, self.thisptr.ksize(), kmer_f, kmer_r)
        return self.thisptr.count_kmers_within_radius(kmer_f, kmer_r, radius,
                                                      max_count, NULL)

    def kmer_degree(self, const char* kmer_s):
        return self.thisptr.kmer_degree(kmer_s)

    def find_radius_for_volume(self, const char *kmer, unsigned int max_count,
                                     unsigned int max_radius):
        cdef HashIntoType kmer_f, kmer_r
        _hash(kmer, self.thisptr.ksize(), kmer_f, kmer_r)
        return self.thisptr.find_radius_for_volume(kmer_f, kmer_r, max_count,
                                                   max_radius)

    def count_kmers_on_radius(self, const char *kmer, unsigned int radius,
                                    unsigned int max_volume=0):
        cdef HashIntoType kmer_f, kmer_r

        _hash(kmer, self.thisptr.ksize(), kmer_f, kmer_r)
        return self.thisptr.count_kmers_on_radius(kmer_f, kmer_r, radius,
                                                  max_volume)

    def add_tag(self, const char *kmer_s):
        cdef HashIntoType kmer = _hash(kmer_s, self.thisptr.ksize())
        self.thisptr.add_tag(kmer)

    def add_stop_tag(self, const char *kmer_s):
        cdef HashIntoType kmer = _hash(kmer_s, self.thisptr.ksize())
        self.thisptr.add_stop_tag(kmer)

    def save_tagset(self, string filename):
        self.thisptr.save_tagset(filename)

    def load_tagset(self, string filename, bool clear_tags=True):
        self.thisptr.load_tagset(filename, clear_tags)

    def consume_fasta_and_tag_with_stoptags(self, const string filename, callback_obj=None):
        cdef unsigned long long n_consumed
        cdef unsigned int total_reads
        cdef CallbackFn _report_fn = NULL

        #try {
        self.thisptr.consume_fasta_and_tag_with_stoptags(filename,
                        total_reads, n_consumed,
                        _report_fn, <void*>callback_obj)
        #} catch (_khmer_signal &e) {
        #    return NULL;
        #}

        return total_reads, n_consumed

    def identify_stoptags_by_position(self, const string seq):
        cdef vector[unsigned int] posns
        self.thisptr.identify_stop_tags_by_position(seq, posns)
        return posns

    def ksize(self):
        return self.thisptr.ksize()

    def hashsizes(self):
        return self.thisptr.get_tablesizes()

    def extract_unique_paths(self, string seq, unsigned int min_length,
                                   float min_unique):
        cdef vector[string] results
        self.thisptr.extract_unique_paths(seq, min_length, min_unique, results)
        return results

    def find_unpart(self, const string filename, bool traverse,
                          bool stop_big_traversals, callback_obj=None):
        cdef unsigned int n_singletons = 0
        cdef CallbackFn _report_fn = NULL

        #try {
        n_singletons = self.thisptr.partition.find_unpart(
                                filename, traverse,
                                stop_big_traversals,
                                _report_fn,<void*>callback_obj)
        #} catch (_khmer_signal &e) {
        #    return NULL;
        #}

        return n_singletons

    def get_median_count(self, const string kmer):
        if len(kmer) < self.thisptr.ksize():
            raise ValueError("k-mer length must be >= the hashtable k-size")

        cdef BoundedCounterType med = 0
        cdef float average = 0, stddev = 0

        self.thisptr.get_median_count(kmer, med, average, stddev)

        return med, average, stddev


_Hashbits = Hashbits


cdef class PrePartitionInfo:
    cdef pre_partition_info *thisptr

    def __dealloc__(self):
        del self.thisptr


cdef PrePartitionInfo ppi_factory(pre_partition_info *o):
    cdef PrePartitionInfo obj = PrePartitionInfo.__new__(PrePartitionInfo)
    obj.thisptr = o
    return obj


cdef class SubsetPartitionInfo:
    cdef SubsetPartition *thisptr

    def __dealloc__(self):
        del self.thisptr


cdef SubsetPartitionInfo spi_factory(SubsetPartition *o):
    cdef SubsetPartitionInfo obj = SubsetPartitionInfo.__new__(SubsetPartitionInfo)
    obj.thisptr = o
    return obj


cdef class ReadParser:
    pass


cdef class KCountingHash:
    cdef CountingHash *thisptr

    # FIXME: number_of_threads should be based in khmer_config!
    def __cinit__(self, WordLength ksize, tablesizes, uint32_t number_of_threads=1):
        try:
            self.thisptr = new CountingHash(ksize,
                    <vector[unsigned long long int]&>tablesizes,
                    number_of_threads)
        except TypeError:
            self.thisptr = new CountingHash(ksize,
                    <HashIntoType>tablesizes,
                    number_of_threads)

    def __dealloc__(self):
        del self.thisptr

    def consume(self, const string& long_str):
        return self.thisptr.consume_string(long_str)

    def get(self, val):
        if isinstance(val, int) or isinstance(val, long):
            return self.thisptr.get_count(<HashIntoType>val);
        else:
            return self.thisptr.get_count(<char *>val);

    def get_median_count(self, const string kmer):
        if len(kmer) < self.thisptr.ksize():
            raise ValueError("k-mer length must be >= the hashtable k-size")

        cdef BoundedCounterType med = 0
        cdef float average = 0, stddev = 0

        self.thisptr.get_median_count(kmer, med, average, stddev)

        return med, average, stddev

    def get_kadian_count(self, const string kmer, unsigned int nk=1):
        if len(kmer) < self.thisptr.ksize():
            raise ValueError("k-mer length must be >= the hashtable k-size")

        cdef BoundedCounterType kad = 0
        self.thisptr.get_kadian_count(kmer, kad, nk)

        return kad

    def consume_fasta(self, const string filename, callback_obj=None):
        cdef unsigned long long n_consumed = 0
        cdef unsigned int total_reads = 0
        cdef CallbackFn _report_fn = NULL

#        try {
        self.thisptr.consume_fasta(filename, total_reads, n_consumed,
                                   _report_fn, <void*>callback_obj)
#        } catch (_khmer_signal &e) {
#            return NULL;
#        }

        return total_reads, n_consumed

    def save(self, const string filename):
        self.thisptr.save(filename)

    def load(self, const string filename):
        self.thisptr.load(filename)

    def abundance_distribution(self, string filename, Hashbits tracking):
        cdef HashIntoType *dist
        dist = self.thisptr.abundance_distribution(filename, tracking.thisptr)
        return [dist[i] for i in range(0, MAX_BIGCOUNT_C)]

    def trim_on_abundance(self, string seq, BoundedCounterType min_count):
        cdef unsigned int trim_at
        trim_at = self.thisptr.trim_on_abundance(seq, min_count)
        return seq.substr(0, trim_at), trim_at

    def set_use_bigcount(self, bool choice):
        self.thisptr.set_use_bigcount(choice)

    def count(self, const char* kmer):
        """Count the given kmer"""
        if len(kmer) != self.thisptr.ksize():
            raise ValueError("k-mer length must be the same as the hashtable k-size")
        self.thisptr.count(kmer)
        return 1

    def ksize(self):
        return self.thisptr.ksize()

    def hashsizes(self):
        return self.thisptr.get_tablesizes()

    def consume_high_abund_kmers(self, string long_str,
                                       BoundedCounterType min_count):
        if len(long_str) < self.thisptr.ksize():
            raise ValueError("string length must >= the hashtable k-mer size")

        if min_count > MAX_COUNT_C:
            raise ValueError("min count specified is > maximum possible count")

        return self.thisptr.consume_high_abund_kmers(long_str, min_count)

    def fasta_count_kmers_by_position(self, string inputfile,
                                            int max_read_len,
                                            int limit_by=0,
                                            callback_obj=None):
        cdef unsigned long long * counts
        cdef CallbackFn _report_fn = NULL

        counts = self.thisptr.fasta_count_kmers_by_position(
                    inputfile,
                    max_read_len,
                    limit_by,
                    _report_fn, <void*>callback_obj)

        return [counts[i] for i in range(0, max_read_len)]

    def get_max_count(self, string long_str):
        if len(long_str) < self.thisptr.ksize():
            raise ValueError("string length must >= the hashtable k-mer size")
        return self.thisptr.get_max_count(long_str)

    def get_min_count(self, string long_str):
        if len(long_str) < self.thisptr.ksize():
            raise ValueError("string length must >= the hashtable k-mer size")
        return self.thisptr.get_min_count(long_str)

    def n_occupied(self, HashIntoType start=0, HashIntoType stop=0):
        return self.thisptr.n_occupied(start, stop)

    def output_fasta_kmer_pos_freq(self, string infile, string outfile):
        self.thisptr.output_fasta_kmer_pos_freq(infile, outfile)
        return 0


cdef class ReadAligner:
    cdef Aligner *thisptr

    def __cinit__(self, KCountingHash ch, lambdaOne=0.0, lambdaTwo=0.0,
                  unsigned int maxErrorRegion=UINT_MAX):
        self.thisptr= new Aligner(ch.thisptr, lambdaOne, lambdaTwo, maxErrorRegion)

    def __dealloc__(self):
        del self.thisptr

    def align(self, string read):
        cdef CandidateAlignment aln
        cdef string rA

        aln = self.thisptr.alignRead(read)
        rA = aln.getReadAlignment(read)

        return aln.alignment, rA


cdef class KTable:
    cdef CppKTable *thisptr

    def __cinit__(self, long size):
        self.thisptr = new CppKTable(size)

    cdef __update_ktable__(self, CppKTable *new):
        # FIXME: only used by intersect, there must be a better way
        del self.thisptr
        self.thisptr = new

    def __dealloc__(self):
        del self.thisptr

    def __len__(self):
        return self.thisptr.n_entries()

    def __getitem__(self, val):
        return self.get(val)

    def __contains__(self, const char* kmer):
        if self.thisptr.get_count(kmer) > 0:
            return True
        return False

    def forward_hash(self, const char* kmer):
        """Convert string to int"""
        if len(kmer) != self.thisptr.ksize():
            raise ValueError("k-mer length must be the same as the hashtable k-size")
        return _hash(kmer, self.thisptr.ksize())

    def forward_hash_no_rc(self, const char* kmer):
        """Convert string to int, with no reverse complement handling"""
        if len(kmer) != self.thisptr.ksize():
            raise ValueError("k-mer length must be the same as the hashtable k-size")
        return _hash_forward(kmer, self.thisptr.ksize())

    def reverse_hash(self, HashIntoType val):
        """Convert int to string"""
        return _revhash(val, self.thisptr.ksize())

    def count(self, const char* kmer):
        """Count the given kmer"""
        if len(kmer) != self.thisptr.ksize():
            raise ValueError("k-mer length must be the same as the hashtable k-size")
        self.thisptr.count(kmer)
        return 1

    def consume(self, const string s):
        """Count all k-mers in the given string"""
        if len(s) < self.thisptr.ksize():
            raise ValueError("string length must >= the hashtable k-mer size")
        self.thisptr.consume_string(s)

        return len(s) - self.thisptr.ksize() + 1

    def get(self, val):
        "Get the count for the given k-mer"
        if isinstance(val, int) or isinstance(val, long):
            return self.thisptr.get_count(<HashIntoType>val);
        else:
            return self.thisptr.get_count(<char *>val);

    def max_hash(self):
        """Get the maximum hash value"""
        return self.thisptr.max_hash()

    def n_entries(self):
        """Get the number of possible entries"""
        return self.thisptr.n_entries()

    def ksize(self):
        """Get k"""
        return self.thisptr.ksize()

    def set(self, val, ExactCounterType c):
        """Set counter to a value"""
        if isinstance(val, int) or isinstance(val, long):
            self.thisptr.set_count(<HashIntoType>val, c);
        else:
            self.thisptr.set_count(<char *>val, c);

    def update(self, KTable other):
        """Combine another ktable with this one"""
        self.thisptr.update(deref(other.thisptr))

    def intersect(self, KTable other):
        """
        Create another ktable containing the intersection of two ktables:
        where both ktables have an entry, the counts will be summed.
        """
        # FIXME: there must be a better way
        cdef CppKTable *intersection
        intersection = self.thisptr.intersect(deref(other.thisptr))
        ktable = KTable(1)
        ktable.__update_ktable__(intersection)
        return ktable

    def clear(self):
        """Set all entries to 0."""
        self.thisptr.clear()



def new_ktable(size):
    """ Create an empty ktable """
    return KTable(size)


def get_config():
    """Get active khmer configuration object"""
    pass


def new_hashtable(k, size):
    """Create an empty single-table counting hash"""
    return KCountingHash(k, size)


def _new_counting_hash(k, sizes, n_threads=1):
    """Create an empty counting hash"""
    return KCountingHash(k, sizes, n_threads)


def _new_hashbits(k, sizes):
    """Create an empty hashbits table"""
    return Hashbits(k, sizes)


def new_readaligner(KCountingHash ch, lambdaOne=0.0, lambdaTwo=0.0, unsigned int maxErrorRegion=UINT_MAX):
    """Create a read aligner object"""
    return ReadAligner(ch, lambdaOne, lambdaTwo, maxErrorRegion)


def forward_hash(const char* kmer, WordLength size):
    max_size = 2 ** (sizeof(WordLength) * 7)
    if size >= max_size:
        raise ValueError("k-mer size must be < %s" % max_size)
    if len(kmer) < size:
        raise ValueError("k-mer size must be >= %s" % size)
    return _hash(kmer, size)


def forward_hash_no_rc(const char* kmer, WordLength size):
    max_size = 2 ** (sizeof(WordLength) * 7)
    if size >= max_size:
        raise ValueError("k-mer size must be < %s" % max_size)
    if len(kmer) < size:
        raise ValueError("k-mer size must be >= %s" % size)
    return _hash_forward(kmer, size)


def reverse_hash(HashIntoType val, WordLength size):
    max_size = 2 ** (sizeof(WordLength) * 7)
    if size >= max_size:
        raise ValueError("k-mer size must be < %s" % max_size)
    return _revhash(val, size)


def set_reporting_callback(cb):
    pass
