from libcpp.string cimport string
from libcpp.memory cimport make_shared
from libc.stdint cimport uint64_t
from cython.operator cimport dereference as deref

from _oxli cimport _revhash

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

    @property
    def kmer_f(self):
        return deref(self._this).kmer_f

    @property
    def kmer_r(self):
        return deref(self._this).kmer_r

    @property
    def kmer_u(self):
        return deref(self._this).kmer_u

    @staticmethod
    cdef Kmer wrap(CpKmer * cpkmer, WordLength K):
        cdef Kmer kmer = Kmer()
        kmer._this.reset(cpkmer)
        kmer.kmer = _revhash(kmer.kmer_u, K).decode('utf-8')
        return kmer

    @staticmethod
    cdef Kmer create(HashIntoType tag, WordLength K):
        cdef Kmer kmer = Kmer()
        deref(kmer._this).set_from_unique_hash(tag, K)
        kmer.kmer = _revhash(kmer.kmer_u, K).decode('utf-8')
        return kmer
