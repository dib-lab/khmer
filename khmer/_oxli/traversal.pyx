from libcpp.memory cimport make_shared
from cython.operator cimport dereference as deref

from khmer._oxli.oxli_types cimport *
from khmer._oxli.graphs cimport Hashgraph
from khmer._oxli.hashing cimport Kmer


cdef class Traverser:

    def __cinit__(self, Hashgraph graph):
        self._graph_ptr = graph._hg_this
        self.graph = graph
        if type(self) is Traverser:
            self._this = make_shared[CpTraverser](self._graph_ptr.get())

    @property
    def ksize(self):
        return self.graph.ksize()

    cdef list _kmerqueue_to_kmer_list(self, KmerQueue * kmers):
        cdef list result = []
        cdef CpKmer cpkmer
        cdef Kmer kmer
        while(deref(kmers).empty() == 0):
            cpkmer = deref(kmers).front()
            kmer = Kmer.wrap(new CpKmer(cpkmer), deref(self._graph_ptr).ksize())
            result.append(kmer)
            deref(kmers).pop()
        return result

    cdef list _kmerqueue_to_hash_list(self, KmerQueue * kmers):
        cdef list result = []
        cdef CpKmer cpkmer
        while(deref(kmers).empty() == 0):
            cpkmer = deref(kmers).front()
            result.append(cpkmer.kmer_u)
            deref(kmers).pop()
        return result

    cdef list _neighbors(self, CpKmer start, int direction=0):
        cdef KmerQueue kmer_q

        if direction == 1:
            deref(self._this).traverse_right(start, kmer_q)
        elif direction == 2:
            deref(self._this).traverse_left(start, kmer_q)
        else:
            deref(self._this).traverse(start, kmer_q)

        cdef list neighbors = self._kmerqueue_to_kmer_list(&kmer_q)

        cdef Kmer neighbor
        for neighbor in neighbors:
            if neighbor.is_forward != start.is_forward():
                neighbor.reverse_complement()
        return neighbors


    def neighbors(self, str node):
        cdef CpKmer start = self.graph._build_kmer(node)
        cdef Kmer neighbor
        for neighbor in self._neighbors(start):
            yield neighbor

    def right_neighbors(self, str node):
        cdef CpKmer start = self.graph._build_kmer(node)
        cdef Kmer neighbor
        for neighbor in self._neighbors(start, direction=1):
            yield neighbor

    def left_neighbors(self, str node):
        cdef CpKmer start = self.graph._build_kmer(node)
        cdef Kmer neighbor
        for neighbor in self._neighbors(start, direction=2):
            yield neighbor

    def degree(self, str node):
        cdef CpKmer kmer = self.graph._build_kmer(node)
        return deref(self._this).degree(kmer)
    
    def left_degree(self, str node):
        cdef CpKmer kmer = self.graph._build_kmer(node)
        return deref(self._this).degree_left(kmer)

    def right_degree(self, str node):
        cdef CpKmer kmer = self.graph._build_kmer(node)
        return deref(self._this).degree_right(kmer)
