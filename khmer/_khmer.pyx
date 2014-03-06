'''
This file is part of khmer, http://github.com/ged-lab/khmer/, and is
Copyright (C) Michigan State University, 2009-2014. It is licensed under
the three-clause BSD license; see doc/LICENSE.txt.
Contact: khmer-project@idyll.org
'''

from cython.operator cimport dereference as deref
from libc.limits cimport UINT_MAX


cdef class _LabelHash:
    pass


cdef class _Hashbits:
    pass


cdef class ReadParser:
    pass


cdef class KCountingHash:
    cdef CountingHash *thisptr

    def __cinit__(self, WordLength ksize, vector[unsigned long long int]& tablesizes,
                  uint32_t number_of_threads):
        self.thisptr = new CountingHash(ksize, tablesizes, number_of_threads)

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


def new_hashtable():
    """Create an empty single-table counting hash"""
    pass


def _new_counting_hash(k, sizes, n_threads=1):
    """Create an empty counting hash"""
    return KCountingHash(k, sizes, n_threads)


def _new_hashbits():
    """Create an empty hashbits table"""
    pass


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
