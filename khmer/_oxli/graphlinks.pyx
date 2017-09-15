cimport cython

from cython.operator cimport dereference as deref

from libcpp.memory cimport shared_ptr, make_shared
from libcpp.string cimport string
from libcpp.vector cimport vector

from khmer._oxli.oxli_types cimport *
from khmer._oxli.hashing cimport CpKmer, CpKmerIterator, _revcomp
from khmer._oxli.sequence cimport Sequence, _object_to_string
from khmer._oxli.sequence import Sequence
from khmer._oxli.traversal cimport CpTraverser
from khmer._oxli.utils cimport _bstring


@cython.freelist(100)
cdef class Junction:

    def __init__(self, HashIntoType u,
                       HashIntoType v):
        self.u = u
        self.v = v

    @staticmethod
    cdef Junction _new(HashIntoType u, 
                       HashIntoType v):
        cdef Junction junc = Junction.__new__(Junction)
        junc.u = u
        junc.v = v
        return junc


cdef class Link:

    def __init__(self, list junctions, bool forward):
        self.junctions = junctions
        self.forward = forward

    @staticmethod
    cdef Link _new(list junctions, bool forward):
        cdef Link link = Link.__new__(Link)
        link.junctions = junctions
        link.forward = forward
        return link


cdef class GraphLinker:

    def __cinit__(self, Hashgraph graph not None):
        self.graph = graph
        self._graph = graph._hg_this
        self.K = graph.ksize()

        self.links = {}

    def check_sequence_length(self, object sequence):
        if len(sequence) < self.K:
            raise ValueError("sequence length ({}) must >= the hashtable "
                             "k-mer size ({})".format(len(sequence),
                                                      self.K))

    cdef tuple _get_junctions(self, string sequence):
        cdef shared_ptr[CpTraverser] _traverser = make_shared[CpTraverser](self._graph.get())
        cdef CpKmerIterator * _it = new CpKmerIterator(sequence.c_str(),
                                                       self.graph.ksize())
        cdef CpKmer u = deref(_it).next()
        cdef CpKmer start = u
        if deref(_it).done():
            return [], []

        cdef list fw_choices = []
        cdef list rc_choices = []
        cdef CpKmer v = deref(_it).next()
        #cdef uint64_t kmer_start_idx = 0
        #cdef uint64_t kmer_end_idx = self.K - 1

        while not deref(_it).done():
            if deref(_traverser).degree_right(u) > 1:
                fw_choices.append(Junction._new(<HashIntoType>u, <HashIntoType>v))
            if deref(_traverser).degree_left(v) > 1:
                rc_choices.append(Junction._new(<HashIntoType>v, <HashIntoType>u))

            u = v
            v = deref(_it).next()
            #kmer_start_idx += 1
            #kmer_end_idx += 1
        rc_choices.reverse()
        return fw_choices, rc_choices

    def get_junctions(self, object sequence):
        self.check_sequence_length(sequence)
        return self._get_junctions(_object_to_string(sequence))

    cdef int _add_link(self, string sequence, unsigned int min_link_size=1) except -1:

        cdef Link fw_link, rc_link
        fw_choices, rc_choices = self._get_junctions(sequence)

        if len(fw_choices) >= min_link_size:
            fw_link = Link._new(fw_choices, True)
            self.links[fw_choices[0].u] = fw_link
        if len(rc_choices) >= min_link_size:
            rc_link = Link._new(rc_choices, False)
            self.links[rc_choices[0].u] = rc_link

        return len(fw_choices) + len(rc_choices)

    def add_link(self, object sequence, unsigned int min_link_size=1):
        if min_link_size < 1:
            raise ValueError("Minimum link size must be at least 1.")
        self.check_sequence_length(sequence)

        return self._add_links(_object_to_string(sequence), min_link_size)

