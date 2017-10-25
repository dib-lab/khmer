# -*- coding: UTF-8 -*-
from khmer._oxli.hashing cimport Kmer, _hash_forward
from khmer._oxli.oxli_types cimport *
from khmer._oxli.utils cimport is_str, is_num, _bstring


cdef class HashSet:
    def __cinit__(self, ksize, hashes=[]):
        self.ksize = ksize
        if hashes:
            self.update(hashes)

    def update(self, hashes):
        for h in hashes:
            self.hs.insert(h)

    cpdef add(self, HashIntoType h):
        self.hs.insert(h)

    cpdef remove(self, HashIntoType h):
        if not self.hs.erase(h):
            raise ValueError('hash not in HashSet; cannot remove')

    def __len__(self):
        return self.hs.size()

    def __contains__(self, object kmer):
        cdef HashIntoType h

        if is_num(kmer):
            h = kmer
        elif isinstance(kmer, Kmer):
            h = kmer.kmer_u
        elif is_str(kmer):
            h = _hash_forward(_bstring(kmer), self.ksize)

        return self.hs.find(h) != self.hs.end()

    def __add__(self, HashSet other):
        if self.ksize != other.ksize:
            raise ValueError('cannot concatenate HashSets with different ksize')
        no = HashSet(self.ksize)
        no += self
        no += other

        return no

    def __iadd__(self, HashSet other):
        cdef HashIntoType h
        if self.ksize != other.ksize:
            raise ValueError('cannot concatenate HashSets with different ksize')
        for h in other.hs:
            self.hs.insert(h)
        return self


    def __iter__(self):
        for h in self.hs:
            yield h
