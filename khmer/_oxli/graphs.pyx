# cython: c_string_type=unicode, c_string_encoding=utf8
from cython.operator cimport dereference as deref

from libcpp.memory cimport unique_ptr

from .utils cimport _bstring
from .._khmer import Countgraph, Nodegraph, GraphLabels


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


cdef class QFCounttable_:
    cdef unique_ptr[CpQFCounttable] c_table
    def __cinit__(self, int k, int size):
        self.c_table.reset(new CpQFCounttable(k, size))

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
        if isinstance(kmer, (unicode, str)):
            data = _bstring(kmer)
            return deref(self.c_table).add(data)
        # assume kmer is an integer representing the hash value
        else:
            return deref(self.c_table).add(kmer)

    def hash(self, kmer):
        """"Returns the hash of this k-mer.

        For Nodetables and Counttables, this function will fail if the
        supplied k-mer contains non-ACGT characters.
        """
        data = _bstring(kmer)
        return deref(self.c_table).hash_dna(data)

    def get(self, kmer):
        """"Retrieve the count for the given k-mer.

        `kmer` can be either a string or an integer representing the hashed
        value of the kmer.

        For Nodetables and Counttables, this function will fail if the
        supplied k-mer contains non-ACGT characters.
        """
        if isinstance(kmer, (unicode, str)):
            data = _bstring(kmer)
            return deref(self.c_table).get_count(data)
        # assume kmer is an integer representing the hash value
        else:
            return deref(self.c_table).get_count(kmer)
