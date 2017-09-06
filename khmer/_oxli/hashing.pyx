# -*- coding: UTF-8 -*-

from libcpp.string cimport string
from libcpp.memory cimport make_shared
from libc.stdint cimport uint64_t
from cython.operator cimport dereference as deref

from khmer._oxli.oxli_types cimport *
from khmer._oxli.utils cimport _bstring, _ustring


cdef class Kmer:

    def __cinit__(self, str kmer=None):
        self.kmer = kmer
        if self.kmer is not None:
            self._this.reset(new CpKmer(kmer.encode('utf-8'), len(kmer)))
        else:
            self._this.reset(new CpKmer())

    def __len__(self):
        return len(self.kmer)

    def __str__(self):
        return self.kmer

    def __hash__(self):
        return self.kmer_u

    def __repr__(self):
        return self.kmer

    @property
    def kmer_f(self):
        return deref(self._this).kmer_f

    @property
    def kmer_r(self):
        return deref(self._this).kmer_r

    @property
    def kmer_u(self):
        return deref(self._this).kmer_u

    def reverse_complement(self):
        cdef HashIntoType tmp = deref(self._this).kmer_f
        deref(self._this).kmer_f = deref(self._this).kmer_r
        deref(self._this).kmer_r = tmp
        self.kmer = _revcomp(self.kmer.encode('utf-8'))

    @property
    def is_forward(self):
        return deref(self._this).is_forward()

    @staticmethod
    cdef Kmer wrap(CpKmer * cpkmer, WordLength K):
        cdef Kmer kmer = Kmer()
        kmer._this.reset(cpkmer)
        kmer.kmer = _revhash(kmer.kmer_u, K)
        return kmer

    @staticmethod
    def create(HashIntoType tag, WordLength K):
        cdef Kmer kmer = Kmer()
        deref(kmer._this).set_from_unique_hash(tag, K)
        kmer.kmer = _revhash(kmer.kmer_u, K)
        return kmer


cpdef HashIntoType forward_hash(str kmer, unsigned int K):
    '''Run the 2-bit hash algorithm on the given K-mer.'''

    if K > 32:
        raise ValueError("k-mer size must be <= 32")
    if len(kmer) != K:
        raise ValueError("k-mer length must equal K")

    return _hash(_bstring(kmer), K)


cpdef HashIntoType forward_hash_no_rc(str kmer, WordLength K):
    '''Run the 2-bit hash function in only the given
    sequence orientation.'''

    if K > 32:
        raise ValueError("k-mer size must be <= 32")
    if len(kmer) != K:
        raise ValueError("k-mer length must equal K")

    return _hash_forward(_bstring(kmer), K)


cpdef str reverse_hash(object h, int K):
    if K > 32:
        raise ValueError("k-mer size must be <= 32")
    
    cdef HashIntoType _h = <HashIntoType>h
    return _revhash(_h, K)


cpdef str reverse_complement(str sequence):
    cdef string s = _revcomp(_bstring(sequence))
    return s


cpdef hash_murmur3(str s):
    cdef HashIntoType h = _hash_murmur(_bstring(s), len(s))
    return h


cpdef hash_no_rc_murmur3(str s):
    cdef HashIntoType h = _hash_murmur_forward(_bstring(s), len(s))
    return h
