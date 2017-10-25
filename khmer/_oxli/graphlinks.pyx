from cython.operator cimport dereference as deref
from libcpp.memory cimport make_shared

from khmer._oxli.utils cimport _bstring
from khmer._oxli.sequence cimport Alphabets


cdef class CompactEdge:

    @staticmethod
    cdef CompactEdge _wrap(CpCompactEdge* _edge):
        cdef CompactEdge edge = CompactEdge()
        edge._ce_this = _edge
        return edge

    def tags(self):
        cdef HashIntoType tag
        for tag in deref(self._ce_this).tags:
            yield tag

    @property
    def edge_type(self):
        cdef compact_edge_meta_t meta = deref(self._ce_this).meta
        if meta == IS_FULL_EDGE:
            return 'FULL'
        elif meta == IS_IN_TIP:
            return 'IN_TIP'
        elif meta == IS_OUT_TIP:
            return 'OUT_TIP'
        elif meta == IS_ISLAND:
            return 'ISLAND'
        elif meta == IS_TRIVIAL:
            return 'TRIVIAL'
        else:
            raise ValueError('Malformed edge metadata')

    @property
    def sequence(self):
        return deref(self._ce_this).sequence

    def in_node(self):
        if deref(self._ce_this).meta == IS_IN_TIP or \
           deref(self._ce_this).meta == IS_ISLAND:
            return None
        return deref(self._ce_this).in_hash

    def out_node(self):
        if deref(self._ce_this).meta == IS_OUT_TIP or \
           deref(self._ce_this).meta == IS_ISLAND:
            return None
        return deref(self._ce_this).out_hash



cdef class CompactNode:

    def __cinit__(self):
        self.kmer = None

    @staticmethod
    cdef CompactNode _wrap(CpCompactNode* _node):
        cdef CompactNode node = CompactNode()
        node._cn_this = _node
        return node

    @property
    def count(self):
        return deref(self._cn_this).count

    @property
    def node_id(self):
        return deref(self._cn_this).node_id

    @property
    def out_degree(self):
        return deref(self._cn_this).out_degree()

    @property
    def in_degree(self):
        return deref(self._cn_this).in_degree()

    def node_kmer(self, WordLength K):
        if self.kmer is None:
            self.kmer = Kmer.wrap(&deref(self._cn_this).kmer, K)
        return self.kmer

    def out_edges(self):
        cdef string bases = Alphabets._get('DNA_SIMPLE')
        cdef char base
        cdef CpCompactEdge * edge
        for base in bases:
            edge = deref(self._cn_this).get_out_edge(base)
            if edge != NULL:
                return base, CompactEdge._wrap(edge)

    def in_edges(self):
        cdef string bases = Alphabets._get('DNA_SIMPLE')
        cdef char base
        cdef CpCompactEdge * edge
        for base in bases:
            edge = deref(self._cn_this).get_in_edge(base)
            if edge != NULL:
                return base, CompactEdge._wrap(edge)


cdef class StreamingCompactor:

    def __cinit__(self, Hashgraph graph):
        self._graph = graph._hg_this
        
        if type(self) is StreamingCompactor:
            self._sc_this = make_shared[CpStreamingCompactor](self._graph)

    def update(self, str sequence):
        cdef string _sequence = _bstring(sequence)
        deref(self._sc_this).update_compact_dbg(_sequence)

    def sequence_nodes(self, str sequence):
        cdef string _sequence = _bstring(sequence)
        cdef vector[uint64_t] node_ids = deref(self._sc_this).get_compact_node_ids(_sequence)
        cdef HashIntoType node_id
        cdef CpCompactNode* node
        for node_id in node_ids:
            node = deref(self._sc_this).get_compact_node_by_id(node_id)
            if node != NULL:
                yield CompactNode._wrap(node)

    def report(self):
        deref(self._sc_this).report()

    @property
    def n_nodes(self):
        return deref(self._sc_this).n_nodes()

    @property
    def n_edges(self):
        return deref(self._sc_this).n_edges()
