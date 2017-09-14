cimport cython

from cython.operator cimport dereference as deref

from libcpp.memory cimport shared_ptr, make_shared
from libcpp.string cimport string
from libcpp.vector cimport vector

from khmer._oxli.oxli_types cimport *
from khmer._oxli.hashing cimport CpKmer, CpKmerIterator
from khmer._oxli.traversal cimport CpTraverser
from khmer._oxli.utils cimport _bstring


@cython.freelist(100)
cdef class Link:

    def __init__(self, HashIntoType u,
                       HashIntoType v,
                       bool forward,
                       Link parent=None):
        self.u = u
        self.v = v
        self.forward = forward
        self.children = [None, None, None, None]
        self.parent = parent

    @cython.boundscheck(False)
    cpdef bool add_child(self, DBGNucl nuc, Link child_link):
        self.children[nuc] = child_link

    @cython.boundscheck(False)
    cpdef Link get_child(self, DBGNucl nuc):
        return self.children[nuc]

    @staticmethod
    cdef Link _new(HashIntoType u, 
                   HashIntoType v, 
                   bool forward,
                   Link parent=None):
        cdef Link link = Link.__new__(Link)
        link.u = u
        link.v = v
        link.forward = forward
        link.parent = parent
        return link


cdef class LinkPath:

    def __cinit__(self):
        self.path = []


cdef class GraphLinker:

    def __cinit__(self, Hashgraph graph not None):
        self.graph = graph
        self._graph = graph._hg_this
        self.K = graph.ksize()

        # mapping from high degree flanking nodes to link paths
        self.links = {}

    cdef list _get_junction_choices(self, string sequence):
        cdef shared_ptr[CpTraverser] _traverser = make_shared[CpTraverser](self._graph.get())
        cdef CpKmerIterator * _it = new CpKmerIterator(sequence.c_str(),
                                                       self.graph.ksize())
        cdef CpKmer u = deref(_it).next()
        if deref(_it).done():
            return []

        cdef list choices = []
        cdef CpKmer v = deref(_it).next()
        cdef uint64_t kmer_start_idx = 0
        cdef uint64_t kmer_end_idx = self.K - 1

        while not deref(_it).done():
            if deref(_traverser).degree_left(v) > 1:

                choices.append(Link._new(<HashIntoType>u, <HashIntoType>v))

            u = v
            v = deref(_it).next()
            kmer_start_idx += 1
            kmer_end_idx += 1
        
        return choices

    def get_junction_choices(self, str sequence):
        cdef list choices = self._get_junction_choices(_bstring(sequence))
        return choices

