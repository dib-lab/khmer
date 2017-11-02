cimport cython
from libcpp.memory cimport shared_ptr
from libcpp.list cimport list as stdlist
from libcpp.pair cimport pair
from libcpp.unordered_set cimport unordered_set as uset
from libcpp.unordered_map cimport unordered_map as umap
from libcpp.vector cimport vector
from libc.stdint cimport uint8_t, uint32_t, uint64_t

from khmer._oxli.oxli_types cimport *
from khmer._oxli.hashing cimport CpKmer, Kmer
from khmer._oxli.graphs cimport CpHashgraph, Hashgraph, Nodegraph, Countgraph


cdef extern from "oxli/links.hh" namespace "oxli" nogil:

    ctypedef pair[HashIntoType, uint64_t] HashIDPair
    ctypedef uset[HashIntoType] UHashSet
    ctypedef vector[HashIntoType] HashVector
    ctypedef umap[HashIntoType, uint64_t] HashIDMap

    cdef cppclass CpCompactEdge "oxli::CompactEdge"
    cdef cppclass CpCompactNode "oxli::CompactNode":
        CpKmer kmer
        uint32_t count
        const uint64_t node_id
        string sequence

        CpCompactEdge* in_edges[4]
        CpCompactEdge* out_edges[4]

        CpCompactNode(CpKmer)
        void add_in_edge(const char, CpCompactEdge*)
        bool delete_in_edge(CpCompactEdge*)
        CpCompactEdge* get_in_edge(const char)
        void add_out_edge(const char, CpCompactEdge*)
        bool delete_out_edge(CpCompactEdge*)
        CpCompactEdge* get_out_edge(const char)

        uint8_t degree()
        uint8_t out_degree()
        uint8_t in_degree()


    ctypedef enum compact_edge_meta_t:
        IS_FULL_EDGE
        IS_IN_TIP
        IS_OUT_TIP
        IS_ISLAND
        IS_TRIVIAL

    cdef cppclass CpCompactEdge "oxli::CompactEdge":
        CpKmer in_node
        CpKmer out_node
        UHashSet tags
        compact_edge_meta_t meta
        string sequence

        CpCompactEdge(HashIntoType, HashIntoType)
        CpComapctEdge(HashIntoType, HashIntoType, compact_edge_meta_t)

        void add_tags(UHashSet&)
        string tag_viz(WordLength)
        float tag_density()

    ctypedef vector[CpCompactNode] CompactNodeVector
    ctypedef umap[HashIntoType, CpCompactEdge*] TagEdgeMap
    ctypedef pair[HashIntoType, CpCompactEdge*] TagEdgePair
    ctypedef set[TagEdgePair] TagEdgePairSet

    cdef cppclass CpStreamingCompactor "oxli::StreamingCompactor":
        shared_ptr[CpHashgraph] graph

        CpStreamingCompactor(shared_ptr[CpHashgraph])
        WordLength ksize()
        void report()
        uint64_t n_nodes()
        uint64_t n_edges()

        CpCompactNode* get_compact_node_by_kmer(HashIntoType)
        CpCompactNode* get_compact_node_by_id(uint64_t)
        CpCompactNode* fetch_or_new_compact_node(CpKmer hdn)
        vector[CpCompactNode*] get_compact_nodes(const string&)
        
        CpCompactEdge* get_compact_edge(uint64_t)
        CpCompactEdge* get_tag_edge_pair(uint64_t, TagEdgePair&)
        CpCompactEdge* get_compact_edge(UHashSet&)

        uint64_t update_compact_dbg(const string&)
        uint64_t consume_sequence(const string&)
        uint64_t consume_sequence_and_update(const string&)


cdef class CompactNode:
    cdef CpCompactNode* _cn_this
    cdef public Kmer kmer

    @staticmethod
    cdef CompactNode _wrap(CpCompactNode*)


cdef class CompactEdge:
    cdef CpCompactEdge* _ce_this

    @staticmethod
    cdef CompactEdge _wrap(CpCompactEdge*)


cdef class StreamingCompactor:

    cdef shared_ptr[CpHashgraph] _graph
    cdef shared_ptr[CpStreamingCompactor] _sc_this

