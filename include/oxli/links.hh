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
#include <functional>
#include <memory>
#include <list>
#include <iostream>
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

namespace oxli {

using std::make_shared;
using std::shared_ptr;

typedef std::pair<HashIntoType, uint64_t> HashIDPair;
typedef std::unordered_set<HashIntoType> UHashSet;
typedef std::vector<HashIntoType> HashVector;
typedef std::unordered_map<HashIntoType, uint64_t> HashIDMap;

class CompactEdge;
class CompactNode {
public:
    Kmer kmer;
    uint32_t count;
    const uint64_t node_id;
    CompactEdge* in_edges[4] = {nullptr, nullptr, nullptr, nullptr};
    CompactEdge* out_edges[4] = {nullptr, nullptr, nullptr, nullptr};

    CompactNode(Kmer kmer, uint64_t node_id) : 
        kmer(kmer), count(0), node_id(node_id) {}

    friend bool operator== (const CompactNode& lhs, const CompactNode& rhs) {
        return lhs.node_id == rhs.node_id;
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
                   << " Count=" << node.count << " in_degree=" 
                   << std::to_string(node.in_degree())
                   << " out_degree=" << std::to_string(node.out_degree()) << ">";
            return stream;
    }
};


enum compact_edge_meta_t {
    IS_FULL_EDGE,
    IS_IN_TIP,
    IS_OUT_TIP,
    IS_ISLAND,
    IS_TRIVIAL
};


class CompactEdge {
public:
    Kmer in_node; // left and right HDNs
    Kmer out_node;
    UHashSet tags;
    compact_edge_meta_t meta;
    std::string sequence;

    CompactEdge(Kmer in_node, Kmer out_node) : 
        in_node(in_node), out_node(out_node), meta(IS_FULL_EDGE) {}
    
    CompactEdge(Kmer in_node, Kmer out_node, compact_edge_meta_t meta) :
        in_node(in_node), out_node(out_node), meta(meta) {}

    void add_tags(UHashSet& new_tags) {
        for (auto tag: new_tags) {
            tags.insert(tag);
        }
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

};


typedef std::vector<CompactNode> CompactNodeVector;
typedef std::vector<CompactEdge> CompactEdgeVector;
typedef std::unordered_map<HashIntoType, CompactEdge*> TagEdgeMap;
typedef std::pair<HashIntoType, CompactEdge*> TagEdgePair;
typedef std::set<TagEdgePair> TagEdgePairSet;
typedef std::set<CompactEdge*> CompactEdgeSet;


class StreamingCompactor
{
    friend class CompactEdge;
    friend class CompactNode;

protected:
    // map from HDN hashes to CompactNode IDs
    HashIDMap hdn_ids;
    // linear storage for CompactNodes
    CompactNodeVector compact_nodes;
    // map from tags to CompactEdges
    TagEdgeMap tags_to_edges;

    uint64_t n_sequences_added;
    uint64_t n_compact_edges;
    uint64_t n_compact_nodes;
    uint32_t tag_density;

    // protected linear creation of CompactNode
    // they should never be deleted, so this is straightforward
    CompactNode* new_compact_node(Kmer hdn) {
        CompactNode * v = get_compact_node_by_kmer(hdn);
        if (v == nullptr) {
            compact_nodes.emplace_back(hdn, n_compact_nodes);
            n_compact_nodes++;
            v = &(compact_nodes.back());
            hdn_ids[hdn] = v->node_id;
        }
        return v;
    }

    CompactEdge* new_compact_edge(Kmer left, Kmer right,
                                   compact_edge_meta_t edge_meta,
                                   std::string edge_sequence) {
         CompactEdge* edge = new CompactEdge(left, right, edge_meta);
         pdebug("new compact edge: left=" << left << " right=" << right
                << " sequence=" << edge_sequence);
         edge->sequence = edge_sequence;
         n_compact_edges++;
         return edge;
    }

    void delete_compact_edge(CompactEdge * edge) {
        pdebug("attempt edge delete @" << edge);
        if (edge != nullptr) {
            pdebug("edge not null, proceeding");
            for (auto tag: edge->tags) {
                tags_to_edges.erase(tag);
            }
            delete edge;
            n_compact_edges--;
        }
    }

    void delete_compact_edge(UHashSet& tags) {
        CompactEdge* edge = get_compact_edge(tags);
        delete_compact_edge(edge);
    }

    void delete_compact_edge(HashIntoType tag) {
        CompactEdge* edge = get_compact_edge(tag);
        delete_compact_edge(edge);
    }

public:

    shared_ptr<Hashgraph> graph;
    
    StreamingCompactor(shared_ptr<Hashgraph> graph) :
        graph(graph), n_sequences_added(0), 
        n_compact_edges(0), n_compact_nodes(0)
    {
        tag_density = DEFAULT_TAG_DENSITY;

    }

    WordLength ksize() const {
        return graph->ksize();
    }

    uint64_t n_nodes() const {
        return compact_nodes.size();
    }

    uint64_t n_edges() const {
        return n_compact_edges;
    }

    void report() const {
        std::cout << "StreamingCompactor(@" << this << " with "
            << "Hashgraph @" << graph.get() << ")" << std::endl;
        std::cout << "  * " << hdn_ids.size() << " cDBG nodes (HDNs)" << std::endl;
        std::cout << "  * " << n_compact_edges << " cDBG edges" << std::endl;
        std::cout << "  * " << n_sequences_added << " sequences added" << std::endl;
    }

    CompactNode* get_compact_node_by_kmer(HashIntoType hdn)  {
        auto search = hdn_ids.find(hdn);
        if (search != hdn_ids.end()) {
            uint64_t ID = search->second;
            return &(compact_nodes[ID]);
        }
        return nullptr;
    }

    CompactNode* get_compact_node_by_id(uint64_t id) {
        if (id < compact_nodes.size()) {
            return nullptr;
        }
        return &(compact_nodes[id]);
    }

    CompactNode* fetch_or_new_compact_node(Kmer hdn) {
        CompactNode* v = get_compact_node_by_kmer(hdn);
        if (v != nullptr) {
            v->count += 1;
        } else {
            v = new_compact_node(hdn);
            v->count = 1;
        }
        return v;
    }

    std::vector<uint64_t> get_compact_node_ids(const std::string& sequence)  {
        KmerIterator kmers(sequence.c_str(), graph->ksize());
        std::vector<uint64_t> ids;

        CompactNode* node;

        while(!kmers.done()) {
            Kmer kmer = kmers.next();

            node = get_compact_node_by_kmer(kmer);
            if (node != nullptr) {
                ids.push_back(node->node_id);
            }
        }

        return ids;
    }

    CompactEdge* get_compact_edge(uint64_t tag) {
        auto search = tags_to_edges.find(tag);
        if (search != tags_to_edges.end()) {
            return search->second;
        }
        return nullptr;
    }

    bool get_tag_edge_pair(uint64_t tag, TagEdgePair& pair) {
        auto search = tags_to_edges.find(tag);
        if (search != tags_to_edges.end()) {
            pair = *search;
            return true;
        } else {
            return false;
        }
    }

    CompactEdge* get_compact_edge(UHashSet& tags) {
        CompactEdge * edge = nullptr;
        for (auto tag: tags) {
            edge = get_compact_edge(tag);
            if (edge != nullptr) {
                break;
            }
        }
        return edge;
    }

    KmerFilter get_tag_stopper(TagEdgePair& te_pair,
                               bool& found_tag) {
        KmerFilter stopper = [&] (const Kmer& node) {
            found_tag = get_tag_edge_pair(node, te_pair);
            return found_tag;
        };

        return stopper;
    }

    compact_edge_meta_t deduce_edge_meta(CompactNode* in, CompactNode* out) {
        compact_edge_meta_t edge_meta;
        if (in == nullptr && out == nullptr) {
            edge_meta = IS_ISLAND;
        } else if (out == nullptr) {
            edge_meta = IS_OUT_TIP;
        } else if (in == nullptr) {
            edge_meta = IS_IN_TIP;
        } else {
            edge_meta = IS_FULL_EDGE;
        }
        return edge_meta;
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

    /* Update a compact dbg where there are no induced
     * HDNs */
    void update_compact_dbg_linear(std::string& sequence) {

    }

    uint64_t update_compact_dbg(const std::string& sequence) {

        // first gather up all k-mers that could have been disturbed --
        // k-mers in the read, and the neighbors of the flanking nodes
        KmerIterator kmers(sequence.c_str(), ksize());
        KmerQueue disturbed_kmers;
        Kmer kmer = kmers.next();
        kmer.set_forward();
        CompactingAT<TRAVERSAL_LEFT> lcursor(graph.get(), kmer);
        lcursor.neighbors(disturbed_kmers);
        while(!kmers.done()) {
            kmer = kmers.next();
            kmer.set_forward();
            disturbed_kmers.push_back(kmer);
        }
        CompactingAT<TRAVERSAL_RIGHT> rcursor(graph.get(), kmer);
        rcursor.neighbors(disturbed_kmers);

        // find the induced HDNs in the disturbed k-mers
        KmerSet induced_hdns;
        uint64_t n_updates = 0;
        while(!disturbed_kmers.empty()) {
            Kmer kmer = disturbed_kmers.back();
            disturbed_kmers.pop_back();
            if(lcursor.degree(kmer) > 1 || rcursor.degree(kmer) > 1) {
                if (fetch_or_new_compact_node(kmer)->count == 1) {
                    induced_hdns.insert(kmer);
                    n_updates++;
                }
            }
        }

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
            root_kmer.set_forward();

            CompactNode* root_node = get_compact_node_by_kmer(root_kmer);

            // check left (in) edges
            lcursor.neighbors(root_kmer, neighbors);
            while(!neighbors.empty()) {
                Kmer neighbor = neighbors.back();
                neighbors.pop_back();
                lcursor.cursor = neighbor;

                TagEdgePair left_tag_pair;
                bool found_left_tag;

                lcursor.push_filter(get_tag_stopper(left_tag_pair, found_left_tag));
                std::string segment_seq = cassem._assemble_directed(lcursor);

                CompactEdge* segment_edge = nullptr;
                CompactNode* left_node = nullptr;
                /*
                 * Should also keep a set of pair<Kmer,Kmer> to track resolved
                 * segments
                 */
                if (found_left_tag) {
                    // must be an existing segment
                    segment_edge = get_compact_edge(left_tag_pair.first);
                    left_node = get_compact_node_by_kmer(segment_edge->in_node);
                    // check if we need to split the edge
                    if (*get_compact_node_by_kmer(segment_edge->out_node) == *root_node) {
                        // good, we're set
                        continue; // continue through neighbors
                    } else {
                        // need to be careful here, the out node for the 
                        // edge we delete could be linked to another induced
                        // HDN...
                        n_updates++;
                        get_compact_node_by_kmer(segment_edge->out_node)->delete_in_edge(segment_edge);
                        delete_compact_edge(segment_edge);
                        // deleted tags; assemble out to the HDN
                        segment_seq = cassem._assemble_directed(lcursor) +
                                      segment_seq.substr(ksize());
                        if (left_node != nullptr) {
                            // not an IN_TIP
                            left_node->delete_out_edge(segment_edge);
                        }
                    }
                } else {
                    // then left node must have been new, or does not exist
                    left_node = get_compact_node_by_kmer(lcursor.cursor);
                }

                // construct the compact edge
                segment_seq = segment_seq.substr(1);
                compact_edge_meta_t edge_meta = (left_node == nullptr) ?
                                                  IS_IN_TIP : IS_FULL_EDGE;

                n_updates++;
                segment_edge = new_compact_edge(lcursor.cursor, 
                                                root_node->kmer,
                                                edge_meta, 
                                                segment_seq);
                if (IS_FULL_EDGE) {
                    left_node->add_out_edge(segment_seq[ksize()], segment_edge);
                } 

            }

            // now the right neighbors...


        }

    }

};



}


#endif
