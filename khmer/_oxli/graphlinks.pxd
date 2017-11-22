cimport cython
from libcpp.memory cimport shared_ptr
from libcpp.list cimport list as stdlist
from libcpp.pair cimport pair
from libcpp.unordered_set cimport unordered_set as uset
from libcpp.unordered_map cimport unordered_map as umap
from libcpp.vector cimport vector
from libc.stdint cimport uint8_t, uint32_t, uint64_t

from khmer._oxli.oxli_types cimport *
from khmer._oxli.hashing cimport CpKmer, Kmer, CpKmerFactory
from khmer._oxli.graphs cimport CpHashgraph, Hashgraph, Nodegraph, Countgraph


cdef extern from "oxli/links.hh":
    cdef uint64_t NULL_ID

cdef extern from "oxli/links.hh" namespace "oxli" nogil:

    ctypedef uint64_t id_t
    ctypedef pair[HashIntoType, id_t] HashIDPair
    ctypedef uset[HashIntoType] UHashSet
    ctypedef vector[HashIntoType] HashVector
    ctypedef umap[HashIntoType, id_t] HashIDMap

    ctypedef enum compact_edge_meta_t:
        FULL
        TIP
        ISLAND
        TRIVIAL

    cdef const char * edge_meta_repr(compact_edge_meta_t)

    cdef cppclass CpCompactEdge "oxli::CompactEdge":
        const id_t in_node_id
        const id_t out_node_id
        const id_t edge_id
        UHashSet tags
        compact_edge_meta_t meta
        string sequence

        CpCompactEdge(id_t, id_t)
        CpComapctEdge(id_t, id_t, compact_edge_meta_t)

        string rc_sequence()
        void add_tags(UHashSet&)
        string tag_viz(WordLength)
        float tag_density()

    ctypedef pair[HashIntoType, CpCompactEdge*] TagEdgePair
    ctypedef set[TagEdgePair] TagEdgePairSet

    cdef cppclass CpCompactEdgeFactory "oxli::CompactEdgeFactory" (CpKmerFactory):
        CpCompactEdgeFactory(WordLength)

        uint64_t n_edges()
        uint64_t n_updates()

        CpCompactEdge* build_edge(id_t, id_t, compact_edge_meta_t,
                                  string)
        void delete_edge(CpCompactEdge*)
        void delete_edge(UHashSet&)
        void delete_edge(HashIntoType)
        CpCompactEdge* get_edge(HashIntoType)
        bool get_tag_edge_pair(HashIntoType, TagEdgePair&)
        CpCompactEdge* get_edge(UHashSet&)

    cdef cppclass CpCompactNode "oxli::CompactNode":
        CpKmer kmer
        uint32_t count
        const id_t node_id
        string sequence

        CpCompactEdge* in_edges[4]
        CpCompactEdge* out_edges[4]

        CpCompactNode(CpKmer, id_t)
        CpCompactNode(CpKmer, string, id_t)

        void add_in_edge(const char, CpCompactEdge*)
        bool delete_in_edge(CpCompactEdge*)
        CpCompactEdge* get_in_edge(const char)
        void add_out_edge(const char, CpCompactEdge*)
        bool delete_out_edge(CpCompactEdge*)
        CpCompactEdge* get_out_edge(const char)
        bool delete_edge(const char)

        uint8_t degree()
        uint8_t out_degree()
        uint8_t in_degree()

    ctypedef vector[CpCompactNode] CompactNodeVector

    cdef cppclass CpCompactNodeFactory "oxli::CompactNodeFactory" (CpKmerFactory):
        CpCompactNodeFactory(WordLength)
        uint64_t n_nodes()
        uint64_t n_updates()

        CpCompactNode* build_node(CpKmer)
        CpCompactNode* get_node_by_kmer(HashIntoType)
        CpCompactNode* get_node_by_id(id_t)
        CpCompactNode* get_or_build_node(CpKmer)
        vector[CpCompactNode*] get_nodes(const string&)

        void unlink_edge(CpCompactEdge*)

        bool get_pivot_from_left(CpCompactNode*, string&, char&)
        bool add_edge_from_left(CpCompactNode*, CpCompactEdge*)
        bool get_edge_from_left(CpCompactNode*, CpCompactEdge* &, string&)

        bool get_pivot_from_right(CpCompactNode*, string&, char&)
        bool add_edge_from_right(CpCompactNode*, CpCompactEdge*)
        bool get_edge_from_right(CpCompactNode*, CpCompactEdge* &, string&)

    cdef cppclass CpStreamingCompactor "oxli::StreamingCompactor":
        shared_ptr[CpHashgraph] graph

        CpStreamingCompactor(shared_ptr[CpHashgraph])
        void report()
        uint64_t n_nodes()
        uint64_t n_edges()
        uint64_t n_updates()

        CpCompactNode* get_node_by_kmer(HashIntoType)
        CpCompactNode* get_node_by_id(id_t)
        vector[CpCompactNode*] get_nodes(const string&)
        
        CpCompactEdge* get_edge(HashIntoType)
        bool get_tag_edge_pair(id_t, TagEdgePair&)
        CpCompactEdge* get_edge(UHashSet&)

        uint64_t update_compact_dbg(const string&)
        uint64_t consume_sequence(const string&)
        uint64_t consume_sequence_and_update(const string&)

        void write_gml(string)


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

