/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2015-2016, The Regents of the University of California.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Michigan State University nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
LICENSE (END)

Contact: khmer-project@idyll.org
*/
#ifndef LINKS_HH
#define LINKS_HH

#include <algorithm>
#include <cstdint>
#include <functional>
#include <memory>
#include <limits>
#include <list>
#include <iostream>
#include <sstream>
#include <unordered_set>
#include <unordered_map>

#include "oxli.hh"
#include "kmer_hash.hh"
#include "hashtable.hh"
#include "hashgraph.hh"
#include "kmer_filters.hh"
#include "traversal.hh"
#include "assembler.hh"


# define DEBUG_LINKS
# ifdef DEBUG_LINKS
#   define pdebug(x) do { std::cout << std::endl << "@ " << __FILE__ <<\
                          ":" << __FUNCTION__ << ":" <<\
                          __LINE__  << std::endl << x << std::endl;\
                          } while (0)
# else
#   define pdebug(x) do {} while (0)
# endif

#define complement(ch) ((ch) == 'A' ? 'T' : \
                        (ch) == 'T' ? 'A' : \
                        (ch) == 'C' ? 'G' : 'C')

namespace oxli {

#define NULL_ID ULLONG_MAX

using std::make_shared;
using std::shared_ptr;

typedef std::pair<HashIntoType, uint64_t> HashIDPair;
typedef std::unordered_set<HashIntoType> UHashSet;
typedef std::vector<HashIntoType> HashVector;
typedef std::unordered_map<HashIntoType, uint64_t> HashIDMap;


enum compact_edge_meta_t {
    IS_FULL_EDGE,
    IS_TIP,
    IS_ISLAND,
    IS_TRIVIAL
};


inline const char * edge_meta_repr(compact_edge_meta_t meta) {
    switch(meta) {
        case IS_FULL_EDGE:
            return "FULL EDGE";
        case IS_TIP:
            return "TIP";
        case IS_ISLAND:
            return "ISLAND";
        case IS_TRIVIAL:
            return "TRIVIAL";
    }
}


class CompactEdgeFactory;
class CompactEdge {
    friend class CompactEdgeFactory;

public:

    const uint64_t in_node_id; // left and right HDN IDs
    const uint64_t out_node_id;
    compact_edge_meta_t meta;
    std::string sequence;
    UHashSet tags;

    CompactEdge(uint64_t in_node_id, uint64_t out_node_id) : 
        in_node_id(in_node_id), out_node_id(out_node_id), meta(IS_FULL_EDGE) {}
    
    CompactEdge(uint64_t in_node_id, uint64_t out_node_id, compact_edge_meta_t meta) :
        in_node_id(in_node_id), out_node_id(out_node_id), meta(meta) {}

    void add_tags(UHashSet& new_tags) {
        for (auto tag: new_tags) {
            tags.insert(tag);
        }
    }

    std::string rc_sequence() const {
        return _revcomp(sequence);
    }

    float tag_density() const {
        return (float)sequence.length() / (float)tags.size();
    }

    std::string tag_viz(WordLength K) const {
        uint64_t pos;
        std::string ret = "L=" + std::to_string(sequence.length()) + " ";
        const char * _s = sequence.c_str();

        for (pos = 0; pos < sequence.length() - K + 1; pos++) {
            if (set_contains(tags, _hash(_s+pos, K))) {
                ret += ("(" + std::to_string(pos) + ")");
            }
            ret += sequence[pos];
        }
        return ret;
    }

    friend std::ostream& operator<<(std::ostream& stream,
                                     const CompactEdge& edge) {
            stream << "<CompactEdge in_node_id=" << 
                      std::to_string(edge.in_node_id) << 
                      " out_node_id=" << 
                      std::to_string(edge.out_node_id)
                   << " length=" << edge.sequence.length()
                   << " meta=" << edge_meta_repr(edge.meta)
                   << " n_tags=" << edge.tags.size() << ">";
            return stream;
    }

};

typedef std::vector<CompactEdge> CompactEdgeVector;
typedef std::unordered_map<HashIntoType, CompactEdge*> TagEdgeMap;
typedef std::pair<HashIntoType, CompactEdge*> TagEdgePair;
typedef std::set<TagEdgePair> TagEdgePairSet;
typedef std::set<CompactEdge*> CompactEdgeSet;

class CompactEdgeFactory : public KmerFactory {

protected:

    uint64_t n_compact_edges;
    uint32_t tag_density;

    TagEdgeMap tags_to_edges;
    CompactEdgeVector compact_edges;

public:

    CompactEdgeFactory(WordLength K) :
        KmerFactory(K), n_compact_edges(0) {

        tag_density = DEFAULT_TAG_DENSITY;    
    }

    uint64_t n_edges() const {
        return n_compact_edges;
    }

    CompactEdge* build_edge(uint64_t left_id, uint64_t right_id,
                            compact_edge_meta_t edge_meta,
                            std::string edge_sequence) {
        CompactEdge* edge = new CompactEdge(left_id, right_id, edge_meta);
        pdebug("new compact edge: \n left=" << std::to_string(left_id) 
                << std::endl << " right=" << std::to_string(right_id)
                << std::endl << " sequence   =" << edge_sequence
                << std::endl << " rc_sequence=" << _revcomp(edge_sequence)
                << std::endl << " start   =" << edge_sequence.substr(0, _ksize)
                << std::endl << " rc_start=" << _revcomp(edge_sequence.substr(0, _ksize))
                << std::endl << " end    =" 
                << edge_sequence.substr(edge_sequence.length()-_ksize, _ksize)
                << std::endl << " rc_end =" 
                << _revcomp(edge_sequence.substr(edge_sequence.length()-_ksize, _ksize)));

        edge->sequence = edge_sequence;
        n_compact_edges++;
        return edge;
    }

    void delete_edge(CompactEdge * edge) {
        //pdebug("attempt edge delete @" << edge);
        if (edge != nullptr) {
            pdebug("edge not null, proceeding");
            for (auto tag: edge->tags) {
                tags_to_edges.erase(tag);
            }
            delete edge;
            n_compact_edges--;
        }
    }

    void delete_edge(UHashSet& tags) {
        CompactEdge* edge = get_edge(tags);
        delete_edge(edge);
    }

    void delete_edge(HashIntoType tag) {
        CompactEdge* edge = get_edge(tag);
        delete_edge(edge);
    }

    CompactEdge* get_edge(HashIntoType tag) const {
        //pdebug("get compact edge from tag " << tag);
        auto search = tags_to_edges.find(tag);
        if (search != tags_to_edges.end()) {
            return search->second;
        }
        return nullptr;
    }

    bool get_tag_edge_pair(HashIntoType tag, TagEdgePair& pair) const {
        auto search = tags_to_edges.find(tag);
        if (search != tags_to_edges.end()) {
            pair = *search;
            return true;
        } else {
            return false;
        }
    }

    CompactEdge* get_edge(UHashSet& tags) const {
        CompactEdge * edge = nullptr;
        for (auto tag: tags) {
            edge = get_edge(tag);
            if (edge != nullptr) {
                break;
            }
        }
        return edge;
    }

    KmerFilter get_tag_stopper(TagEdgePair& te_pair,
                               bool& found_tag) { const
        KmerFilter stopper = [&] (const Kmer& node) {
            found_tag = get_tag_edge_pair(node, te_pair);
            return found_tag;
        };

        return stopper;
    }
};


class CompactNodeFactory;
class CompactNode {
    friend class CompactNodeFactory;
public:
    Kmer kmer;
    uint32_t count;
    const uint64_t node_id;
    std::string sequence;
    bool direction;

    CompactEdge* in_edges[4] = {nullptr, nullptr, nullptr, nullptr};
    CompactEdge* out_edges[4] = {nullptr, nullptr, nullptr, nullptr};

    CompactNode(Kmer kmer, uint64_t node_id) : 
        kmer(kmer), count(0), node_id(node_id), direction(kmer.is_forward()) {}

    CompactNode(Kmer kmer, std::string sequence, uint64_t node_id) : 
        kmer(kmer), count(0), sequence(sequence), node_id(node_id),
        direction(kmer.is_forward()) {}

    friend bool operator== (const CompactNode& lhs, const CompactNode& rhs) {
        return lhs.node_id == rhs.node_id;
    }

    bool delete_edge(CompactEdge* edge) {
        if (delete_in_edge(edge)) {
            return true;
        }
        if (delete_out_edge(edge)) {
            return true;
        }
        return false;
    }

    bool delete_in_edge(CompactEdge* edge) {
        for (uint8_t i=0; i<4; i++) {
            if (in_edges[i] == edge) {
                in_edges[i] = nullptr;
                return true;
            }
        }
        return false;
    }

    void add_in_edge(const char base, CompactEdge* edge) {
        //pdebug("add in edge to " << *this << ", base=" << base
        //        << ", edge: " << *edge);
        in_edges[twobit_repr(base)] = edge;
    }

    CompactEdge* get_in_edge(const char base) {
        return in_edges[twobit_repr(base)];
    }

    bool delete_out_edge(CompactEdge* edge) {
        for (uint8_t i=0; i<4; i++) {
            if (out_edges[i] == edge) {
                out_edges[i] = nullptr;
                return true;
            }
        }
        return false;
    }

    void add_out_edge(const char base, CompactEdge* edge) {
        //pdebug("add out edge to " << *this << ", base=" << base
        //        << ", edge: " << *edge);
        out_edges[twobit_repr(base)] = edge;
    }

    CompactEdge* get_out_edge(const char base) {
        return out_edges[twobit_repr(base)];
    }

    uint8_t degree() const {
        return out_degree() + in_degree();
    }

    uint8_t out_degree() const {
        uint8_t acc = 0;
        for (auto edge: out_edges) {
            if (edge != nullptr) {
                acc++;
            }
        }
        return acc;
    }

    uint8_t in_degree() const {
        uint8_t acc = 0;
        for (auto edge: in_edges) {
            if (edge != nullptr) {
                acc++;
            }
        }
        return acc;
    }

    friend std::ostream& operator<<(std::ostream& stream,
                                     const CompactNode& node) {
            stream << "<CompactNode ID=" << node.node_id << " Kmer=" << node.kmer.kmer_u
                   << " Sequence=" << node.sequence
                   << " Count=" << node.count << " in_degree=" 
                   << std::to_string(node.in_degree())
                   << " out_degree=" << std::to_string(node.out_degree()) << ">";
            return stream;
    }

    std::string edges_repr() {
        std::ostringstream os;
        os << *this << std::endl << "\tin_edges:" << std::endl;
        for (const char b : "ACGT") {
            CompactEdge* e = get_in_edge(b);
            if (e != nullptr) {
                os << "\t " << b << "=" << *e << std::endl;
            }
        }
        os << "\tout_edges:" << std::endl;
        for (const char b : "ACGT") {
            CompactEdge* e = get_out_edge(b);
            if (e != nullptr) {
                os << "\t -" << b << "=" << *e << std::endl;
            }
        }
        return os.str();
    }
};

typedef std::vector<CompactNode> CompactNodeVector;

class CompactNodeFactory : public KmerFactory {

protected:

    // map from HDN hashes to CompactNode IDs
    HashIDMap kmer_id_map;
    // linear storage for CompactNodes
    CompactNodeVector compact_nodes;
    uint64_t n_compact_nodes;

public:
    CompactNodeFactory(WordLength K) : 
        KmerFactory(K), n_compact_nodes(0) {}

    uint64_t n_nodes() const {
        return n_compact_nodes;
    }
    

    // protected linear creation of CompactNode
    // they should never be deleted, so this is straightforward
    CompactNode* build_node(Kmer hdn) {
        //pdebug("new compact node from " << hdn);
        CompactNode * v = get_node_by_kmer(hdn);
        if (v == nullptr) {
            compact_nodes.emplace_back(hdn, n_compact_nodes);
            n_compact_nodes++;
            v = &(compact_nodes.back());
            v->sequence = _revhash(hdn, _ksize);
            kmer_id_map[hdn] = v->node_id;
            //pdebug("Allocate: " << *v);
        }
        return v;
    }

    CompactNode* get_node_by_kmer(HashIntoType hdn) {
        auto search = kmer_id_map.find(hdn);
        if (search != kmer_id_map.end()) {
            uint64_t ID = search->second;
            return &(compact_nodes[ID]);
        }
        return nullptr;
    }

    CompactNode* get_node_by_id(uint64_t id) {
        if (id >= compact_nodes.size()) {
            return nullptr;
        }
        return &(compact_nodes[id]);
    }

    CompactNode* get_or_build_node(Kmer hdn) {
        CompactNode* v = get_node_by_kmer(hdn);
        if (v != nullptr) {
            v->count += 1;
        } else {
            v = build_node(hdn);
            v->count = 1;
        }
        return v;
    }

    std::vector<CompactNode*> get_nodes(const std::string& sequence) {
        //pdebug("get compact node IDs");
        KmerIterator kmers(sequence.c_str(), _ksize);
        std::vector<CompactNode*> nodes;

        CompactNode* node;

        while(!kmers.done()) {
            Kmer kmer = kmers.next();

            node = get_node_by_kmer(kmer);
            if (node != nullptr) {
                nodes.push_back(node);
            }
        }

        return nodes;
    }

    uint8_t unlink_edge(CompactEdge* edge) {
        pdebug("unlink edge " << *edge);
        CompactNode *left, *right;
        left = get_node_by_id(edge->in_node_id);
        right = get_node_by_id(edge->out_node_id);
        uint8_t n_deletes = 0;
        if (left != nullptr) {
            // be lazy for now and use bidirectional delete
            left->delete_edge(edge);
            n_deletes++;
        }
        if (right != nullptr) {
            right->delete_edge(edge);
            n_deletes++;
        }
        return n_deletes;
    }

    bool get_pivot_from_left(CompactNode* v,
                                   std::string& sequence,
                                   char& pivot_base) const {
        const char * node_kmer = v->sequence.c_str();
        const char * _segment = sequence.c_str();
        pivot_base = _segment[sequence.size()-_ksize];
        if (strncmp(node_kmer, 
                    _segment+(sequence.size())-_ksize+1, 
                    _ksize-1) == 0) {
            // same canonical orientation
            return false;
        } else {
            // must have opposite canonical orientation
            pivot_base = complement(pivot_base);
            return true;
        }
    }

    bool add_edge_from_left(CompactNode* v, CompactEdge* e) const {
        char pivot_base;
        if (!get_pivot_from_left(v, e->sequence, pivot_base)) {
            // same canonical orientation
            pdebug("add in edge " << *e << " to node " << *v << " from " << pivot_base);
            v->add_in_edge(pivot_base, e);
            return false;
        } else {
            // must have opposite canonical orientation
            pdebug("add out edge " << *e << " to node " << *v << " from " << pivot_base);
            v->add_out_edge(pivot_base, e);
            return true;
        }
    }


    bool get_edge_from_left(CompactNode* v,
                            CompactEdge* &result_edge,
                            std::string& sequence) const {
        char pivot_base;
        if (!get_pivot_from_left(v, sequence, pivot_base)) {
            result_edge = v->get_in_edge(pivot_base);
            return false;
        } else {
            result_edge = v->get_out_edge(pivot_base);
            return true;
        }
    }

    bool get_pivot_from_right(CompactNode* v,
                              std::string& sequence,
                              char& pivot_base) const {
        const char * node_kmer = v->sequence.c_str();
        const char * _segment = sequence.c_str();
        pivot_base = _segment[_ksize-1];
        if (strncmp(node_kmer+1, _segment, _ksize-1) == 0) {
            // same canonical orientation
            return false;
        } else {
            // must have opposite canonical orientation
            pivot_base = complement(pivot_base);
            return true;
        }
    }

    bool add_edge_from_right(CompactNode* v, CompactEdge* e) const {
        char pivot_base;
        if (!get_pivot_from_right(v, e->sequence, pivot_base)) {
            pdebug("add out edge " << *e << " to node " << *v << " from " << pivot_base);
            v->add_out_edge(pivot_base, e);
            return false;
        } else {
            pdebug("add in edge " << *e << " to node " << *v << " from " << pivot_base);
            v->add_in_edge(pivot_base, e);
            return true;
        }
    }

    bool get_edge_from_right(CompactNode* v,
                             CompactEdge* &result_edge,
                             std::string& sequence) const {
        char pivot_base;
        if (!get_pivot_from_right(v, sequence, pivot_base)) {
            result_edge = v->get_out_edge(pivot_base);
            return false;
        } else {
            result_edge = v->get_in_edge(pivot_base);
            return true;
        }

    }
};


class StreamingCompactor : public KmerFactory
{

protected:

    // map from tags to CompactEdges
    CompactNodeFactory nodes;
    CompactEdgeFactory edges;

    uint64_t n_sequences_added;

public:

    shared_ptr<Hashgraph> graph;
    
    StreamingCompactor(shared_ptr<Hashgraph> graph) :
        KmerFactory(graph->ksize()),
        nodes(graph->ksize()), edges(graph->ksize()),
        n_sequences_added(0), graph(graph)
    {
    }

    compact_edge_meta_t deduce_meta(CompactNode* in, CompactNode* out) {
        compact_edge_meta_t edge_meta;
        if (in == nullptr && out == nullptr) {
           edge_meta = IS_ISLAND;
        } else if ((out == nullptr) != (in == nullptr))  {
            edge_meta = IS_TIP;
        } else {
            edge_meta = IS_FULL_EDGE;
        }
        return edge_meta;
    }

    uint64_t n_nodes() const {
        return nodes.n_nodes();
    }

    uint64_t n_edges() const {
        return edges.n_edges();
    }

    void report() const {
        std::cout << std::endl << "REPORT: StreamingCompactor(@" << this << " with "
            << "Hashgraph @" << graph.get() << ")" << std::endl;
        std::cout << "  * " << n_nodes() << " cDBG nodes (HDNs)" << std::endl;
        std::cout << "  * " << n_edges() << " cDBG edges" << std::endl;
        std::cout << "  * " << n_sequences_added << " sequences added" << std::endl;
    }


    CompactNode* get_node_by_kmer(Kmer hdn) {
        return nodes.get_node_by_kmer(hdn);
    }

    CompactNode* get_node_by_id(uint64_t id) {
        return nodes.get_node_by_id(id);
    }

    std::vector<CompactNode*> get_nodes(const std::string& sequence) {
        return nodes.get_nodes(sequence);
    }

    CompactEdge* get_edge(HashIntoType tag) const {
        return edges.get_edge(tag);
    }

    bool get_tag_edge_pair(HashIntoType tag, TagEdgePair& pair) const {
        return edges.get_tag_edge_pair(tag, pair);
    }

    CompactEdge* get_edge(UHashSet& tags) const {
        return edges.get_edge(tags);
    }

    uint64_t consume_sequence(const std::string& sequence) {
        uint64_t prev_n_kmers = graph->n_unique_kmers();
        graph->consume_string(sequence);
        return graph->n_unique_kmers() - prev_n_kmers;
    }

    uint64_t consume_sequence_and_update(const std::string& sequence) {
        if (consume_sequence(sequence) > 0) {
            return update_compact_dbg(sequence);
        }
        return 0;
    }

    bool validate_segment(CompactNode* root_node, CompactNode* other_node,
                          CompactEdge* edge, std::string& sequence) {
        bool edge_invalid;
        if (edge->meta == IS_TIP) {
            if (other_node != nullptr) {
                return false;
            }
            if (!((edge->in_node_id == root_node->node_id ||
                   edge->out_node_id == root_node->node_id) &&
                  edge->sequence.length() == sequence.length())) {
                return false;
            }
        } else if (edge->meta == IS_FULL_EDGE) {
            if (other_node == nullptr) {
                return false;
            } else {
                bool nodes_match;
                nodes_match = (edge->in_node_id == root_node->node_id && 
                               edge->out_node_id == other_node->node_id) ||
                              (edge->out_node_id == root_node->node_id &&
                               edge->in_node_id == other_node->node_id);
                if (!nodes_match) {
                    return false;
                }
            }
        }
        return true;
    }

    /* Update a compact dbg where there are no induced
     * HDNs
    uint64_t update_compact_dbg_linear(std::string& sequence) {
        Kmer root_kmer = graph->build_kmer(sequence.substr(0, _ksize));

        CompactingAT<TRAVERSAL_LEFT> lcursor(graph.get(), root_kmer);
        CompactingAT<TRAVERSAL_RIGHT> rcursor(graph.get(), root_kmer);
        CompactingAssembler cassem(graph.get());

        std::string left_seq = cassem._assemble_directed(lcursor);
        std::string right_seq = cassem._assemble_directed(rcursor);
        std::string segment_seq = left_seq + right_seq.substr(_ksize);

        CompactNode *left_node = nullptr, *right_node = nullptr;
        left_node = nodes.get_node_by_kmer(lcursor.cursor);
        right_node = nodes.get_node_by_kmer(rcursor.cursor);

        compact_edge_meta_t edge_meta = deduce_edge_meta(left_node, right_node);
        char in_base = segment_seq[_ksize-1];
        char out_base = segment_seq[segment_seq.length()-_ksize+1];

        switch(edge_meta) {
            case IS_FULL_EDGE:
                // then we merged two tips
                pdebug("merge TIPs");
                break;
            case IS_TIP:
                pdebug("extend TIP");
                break;
            case IS_ISLAND:
                pdebug("created or extended ISLAND");
                break;
        }
    }
    */

    uint64_t update_compact_dbg(const std::string& sequence) {
        pdebug("update cDBG from " << sequence);
        n_sequences_added++;

        // first gather up all k-mers that could have been disturbed --
        // k-mers in the read, and the neighbors of the flanking nodes
        KmerIterator kmers(sequence.c_str(), _ksize);
        KmerQueue disturbed_kmers;
        Kmer kmer = kmers.next();
        CompactingAT<TRAVERSAL_LEFT> lcursor(graph.get(), kmer);
        lcursor.neighbors(disturbed_kmers);
        while(!kmers.done()) {
            kmer = kmers.next();
            disturbed_kmers.push_back(kmer);
        }
        CompactingAT<TRAVERSAL_RIGHT> rcursor(graph.get(), kmer);
        rcursor.neighbors(disturbed_kmers);
        
        pdebug(disturbed_kmers.size() << " k-mers disturbed" << std::endl);
        
                // find the induced HDNs in the disturbed k-mers
        KmerSet induced_hdns;
        uint64_t n_updates = 0;
        while(!disturbed_kmers.empty()) {
            Kmer kmer = disturbed_kmers.back();
            disturbed_kmers.pop_back();
            uint8_t l_degree, r_degree;
            l_degree = lcursor.degree(kmer);
            r_degree = rcursor.degree(kmer);
            if(l_degree > 1 || r_degree > 1) {
                pdebug("found HDN... " << kmer);
                CompactNode* hdn = nodes.get_or_build_node(kmer);
                if (hdn->count == 1) {
                    induced_hdns.insert(kmer);
                    n_updates++;
                } else if (hdn->degree() != (l_degree + r_degree)) {
                    induced_hdns.insert(kmer);
                    n_updates++;
                }
            }
        }
        pdebug(induced_hdns.size() << " induced HDNs");

        /* If there are no induced HDNs, we must have extended
         * a tip or merged two tips into a linear segment */
        // handle_merge()


        /* Update from all induced HDNs
         */
        CompactingAssembler cassem(graph.get());
        KmerQueue neighbors;
        while(!induced_hdns.empty()) {
            Kmer root_kmer = *induced_hdns.begin();
            induced_hdns.erase(root_kmer);

            CompactNode* root_node = nodes.get_node_by_kmer(root_kmer);
            pdebug("searching from induced HDN: " << root_node->edges_repr());

            // check left (in) edges
            lcursor.neighbors(root_kmer, neighbors);
            pdebug("checking " << neighbors.size() << " left neighbors");
            while(!neighbors.empty()) {
                Kmer neighbor = neighbors.back();
                neighbors.pop_back();
                lcursor.cursor = neighbor;

                TagEdgePair tag_pair;
                bool found_tag = false;

                lcursor.push_filter(edges.get_tag_stopper(tag_pair, found_tag));
                std::string segment_seq = cassem._assemble_directed(lcursor);

                // first check for a segment going this direction from root
                CompactEdge* segment_edge = nullptr;
                bool edge_invalid = false;
                bool root_flipped = nodes.get_edge_from_left(root_node, segment_edge, segment_seq);

                CompactNode* left_node = nodes.get_node_by_kmer(lcursor.cursor);

                // validate edge leaving root if it exists
                if (segment_edge != nullptr && 
                    validate_segment(root_node, left_node, segment_edge, segment_seq)) {

                    pdebug("validated " << *root_node << ", " << *left_node 
                            << ", " << *segment_edge << ", " << segment_seq);
                    continue;
                }
                
                /*
                 * Should also keep a set of pair<Kmer,Kmer> to track resolved
                 * segments
                 */

                // if there's an existing edge, check if we need to split it
                if (segment_edge != nullptr) {

                    // check that it isn't a TIP from an existing node
                    // with ROOT being induced; if not, delete its edge
                    // from the out node
                    pdebug("edge does not conform, delete it");
                    n_updates += nodes.unlink_edge(segment_edge);
                    edges.delete_edge(segment_edge);
                    // deleted tags; assemble out to the HDN
                    segment_seq = cassem._assemble_directed(lcursor) +
                                  segment_seq.substr(_ksize);
                }

                // construct the compact edge
                compact_edge_meta_t edge_meta = (left_node == nullptr) ?
                                                  IS_TIP : IS_FULL_EDGE;

                if (edge_meta == IS_FULL_EDGE) {
                    // left side includes HDN, right side does not

                    segment_edge = edges.build_edge(left_node->node_id, 
                                                    root_node->node_id,
                                                    edge_meta, 
                                                    segment_seq.substr(1));
                    nodes.add_edge_from_right(left_node, segment_edge);
                } else {
                    segment_edge = edges.build_edge(NULL_ID, 
                                                    root_node->node_id,
                                                    edge_meta, 
                                                    segment_seq);
                }


                n_updates++;
                nodes.add_edge_from_left(root_node, segment_edge);
            }

            // now the right neighbors...
            rcursor.neighbors(root_kmer, neighbors);
            pdebug("checking " << neighbors.size() << " right neighbors");
            while(!neighbors.empty()) {
                Kmer neighbor = neighbors.back();
                neighbors.pop_back();
                rcursor.cursor = neighbor;

                TagEdgePair tag_pair;
                bool found_tag = false;

                rcursor.push_filter(edges.get_tag_stopper(tag_pair, found_tag));
                std::string segment_seq = cassem._assemble_directed(rcursor);

                // first check for a segment going this direction from root
                CompactEdge* segment_edge = nullptr;
                bool edge_invalid = false;
                bool root_flipped = nodes.get_edge_from_right(root_node, segment_edge, segment_seq);

                CompactNode* right_node = nodes.get_node_by_kmer(rcursor.cursor);

                // validate edge leaving root if it exists
                if (segment_edge != nullptr && 
                    validate_segment(root_node, right_node, segment_edge, segment_seq)) {

                    continue;
                }

                if (segment_edge != nullptr) {
                    pdebug("edge does not conform, delete it");
                    n_updates += nodes.unlink_edge(segment_edge);
                    edges.delete_edge(segment_edge);
                    // deleted tags; assemble out to the HDN
                    segment_seq = segment_seq + 
                                  cassem._assemble_directed(rcursor).substr(_ksize);
                }

                compact_edge_meta_t edge_meta = (right_node == nullptr) ?
                                                  IS_TIP : IS_FULL_EDGE;

                if (edge_meta == IS_FULL_EDGE) {
                    segment_edge = edges.build_edge(root_node->node_id, 
                                                    right_node->node_id,
                                                    edge_meta, 
                                                    segment_seq.substr(0, segment_seq.length()-1));
                    nodes.add_edge_from_left(right_node, segment_edge);
                } else {
                    segment_edge = edges.build_edge(root_node->node_id,
                                                    NULL_ID,
                                                    edge_meta, 
                                                    segment_seq);
                }


                n_updates++;
                nodes.add_edge_from_right(root_node, segment_edge);
            }

        }

        return n_updates;

    } // update_compact_dbg

};



}


#endif
