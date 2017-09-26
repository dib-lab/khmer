# -*- coding: UTF-8 -*-

from libcpp.string cimport string
from libcpp.memory cimport make_shared
from libc.stdint cimport uint64_t
from cython.operator cimport dereference as deref

from khmer._oxli.oxli_types cimport *
from khmer._oxli.utils cimport _bstring

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


cdef class SequenceHasher:

    def __cinit__(self, str sequence, WordLength K):
        self._sequence = _bstring(sequence)
        self.K = K

        if type(self) is SequenceHasher:
            print("__cinit__ SequenceHasher({K})".format(K=K))
            self._hasher = make_shared[CpTwoBitKmerHasher](K)
            self._sh_this = make_shared[CpSequenceHasher](self._sequence,
                                                          self._hasher.get())

    def __iter__(self):
        while not deref(self._sh_this).done():
            yield deref(self._sh_this).next().key()
