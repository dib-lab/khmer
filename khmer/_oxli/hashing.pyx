from libcpp.string cimport string
from libcpp.memory cimport shared_ptr
from libc.stdint cimport uint64_t
from cython.operator cimport dereference as deref

cimport ._oxli as _ox

cdef class Kmer:

    cdef shared_ptr[_ox.Kmer] _this
    cdef readonly str kmer

    def __cinit__(self, kmer=None):
        self.kmer = kmer
        if self.kmer is not None:
            self._this.reset(new _ox.Kmer(kmer, len(kmer)))
        else:
            self._this.reset(new _ox.Kmer())

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
    cdef Kmer create(_ox.Kmer& cpkmer):
        cdef Kmer kmer = Kmer()
        kmer._this.reset(&cpkmer)
        return kmer

    @staticmethod
    cdef Kmer create(_ox.HashIntoType tag, _ox.WordLength K):
        cdef Kmer kmer = Kmer()
        deref(kmer._this).set_from_unique_hash(tag, K)
        return kmer
