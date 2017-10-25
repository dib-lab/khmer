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

    void add_in_edge(const char base, CompactEdge* edge) {
        in_edges[twobit_repr(base)] = edge;
    }

    CompactEdge* get_in_edge(const char base) {
        return in_edges[twobit_repr(base)];
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

    friend std::ostream& operator<< (std::ostream& stream,
                                     const CompactNode& node);
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
    HashIntoType in_hash; // left and right HDNs
    HashIntoType out_hash;
    UHashSet tags;
    compact_edge_meta_t meta;
    bool dirty;
    std::string sequence; // optional

    CompactEdge(HashIntoType in_hash, HashIntoType out_hash) : 
        in_hash(in_hash), out_hash(out_hash), meta(IS_FULL_EDGE), dirty(false) {}
    
    CompactEdge(HashIntoType in_hash, HashIntoType out_hash, compact_edge_meta_t meta) :
        in_hash(in_hash), out_hash(out_hash), meta(meta), dirty(false) {}

    void add_tags(UHashSet& new_tags) {
        for (auto tag: new_tags) {
            tags.insert(tag);
        }
    }
};


typedef std::vector<CompactNode> CompactNodeVector;
typedef std::unordered_map<HashIntoType, CompactEdge*> TagEdgeMap;
typedef std::pair<HashIntoType, CompactEdge*> TagEdgePair;
typedef std::set<TagEdgePair> TagEdgePairSet;


class StreamingCompactor
{
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

    CompactEdge* new_compact_edge(HashIntoType left, HashIntoType right,
                                   compact_edge_meta_t edge_meta,
                                   std::string edge_sequence) {
         CompactEdge* edge = new CompactEdge(left, right, edge_meta);
         edge->sequence = edge_sequence;
         n_compact_edges++;
         return edge;
    }

    void delete_compact_edge(CompactEdge * edge) {
        if (edge != nullptr) {
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

    KmerFilter get_tag_finder(TagEdgePairSet& found_pairs,
                              UHashSet& found_tags,
                              UHashSet& all_unresolved_tags) {

        KmerFilter finder = [&] (const Kmer& node) {

            TagEdgePair edge_pair;
            bool found_pair = get_tag_edge_pair(node, edge_pair);
            if (found_pair) {
                found_pairs.insert(edge_pair);
                found_tags.insert(node);
            }
            if (set_contains(all_unresolved_tags, node)) {
                found_tags.insert(node);
            }
            return false; // never filter node, just gather while we're looking
        };

        return finder;
    }

    void update_compact_dbg(UHashSet& unresolved_tags) {

        HashVector unresolved_q;
        for (auto tag: unresolved_tags) {
            unresolved_q.push_back(tag);
        }

        CompactingAssembler at(graph.get());
        
        while(unresolved_q.size() > 0) {
            HashIntoType seed_tag = unresolved_q.back();
            Kmer seed_kmer = graph->build_kmer(seed_tag);
            unresolved_q.pop_back();

            if (!set_contains(unresolved_tags, seed_tag)) {
                continue;
            }

            // build filters to prepare for traversal...
            // first, tag gatherer
            TagEdgePairSet segment_edges;
            UHashSet segment_tags;
            KmerFilterList filters;
            filters.push_back(get_tag_finder(segment_edges,
                                             segment_tags,
                                             unresolved_tags));
            // ... and a shared list of known k-mers
            shared_ptr<SeenSet> visited = make_shared<SeenSet>();
            filters.push_back(get_visited_filter(visited));
            CompactingAT<TRAVERSAL_LEFT> lcursor(graph.get(), seed_kmer,
                                                 filters);
            CompactingAT<TRAVERSAL_RIGHT> rcursor(graph.get(), seed_kmer,
                                                  filters);

            // assemble both dirs until HDN, collecting tags along the way
            std::string left_seq = at._assemble_directed(lcursor);
            std::string right_seq = at._assemble_directed(rcursor);
            std::string edge_sequence = left_seq + right_seq.substr(graph->ksize());

            // Generate new CompactNodes if needed
            CompactNode *left_node = nullptr, *right_node = nullptr;
            // Check degree and gather neighbors simultaneously
            // can't just check if total degree >2, as node could be:
            //        _(neighbor_1)
            //  (HDN)/
            //       \_(neighbor_2)
            // which we'd still like in the cDBG
            KmerQueue rneighbors;
            if (rcursor.neighbors(rneighbors) > 1 || 
                lcursor.neighbors(rcursor.cursor, rneighbors) > 1) {

                right_node = fetch_or_new_compact_node(rcursor.cursor);
                if (right_node->count == 1) {
                    while(rneighbors.size() > 0) {
                        // if it was a new HDN, we it could be between tag and
                        // nearest existing HDN; so, we have to check from its
                        // neighbors
                        auto rn = rneighbors.back();
                        rneighbors.pop();
                        unresolved_q.push_back(rn);
                        unresolved_tags.insert(rn);
                    }
                }
            }
            // same for other side
            KmerQueue lneighbors;
            if (lcursor.neighbors(lneighbors) > 1 || 
                rcursor.neighbors(lcursor.cursor, lneighbors) > 1) {

                left_node = fetch_or_new_compact_node(lcursor.cursor);
                if (left_node->count == 1) {
                    while(lneighbors.size() > 0) {
                        auto ln = lneighbors.back();
                        lneighbors.pop();
                        unresolved_q.push_back(ln);
                        unresolved_tags.insert(ln);
                    }
                }
            }

            // Determine if we've got a tip, full edge, or island
            compact_edge_meta_t edge_meta;
            if (left_node == nullptr && right_node == nullptr) {
                edge_meta = IS_ISLAND;
            } else if (right_node == nullptr) {
                edge_meta = IS_OUT_TIP;
            } else if (left_node == nullptr) {
                edge_meta = IS_IN_TIP;
            } else {
                edge_meta = IS_FULL_EDGE;
            }

            // first check if we've got a good Edge to work with
            CompactEdge* consensus_edge = nullptr;
            TagEdgePair consensus_tag_edge_pair;
            for (auto tag_edge_pair: segment_edges) {
                CompactEdge* intersect_edge = tag_edge_pair.second;
                if (intersect_edge->meta != edge_meta) {
                    continue;
                }

                if (intersect_edge->in_hash == (HashIntoType)(left_node->kmer) &&
                    intersect_edge->out_hash == (HashIntoType)(right_node->kmer)) {
                    
                    consensus_edge = intersect_edge;
                    consensus_tag_edge_pair = tag_edge_pair;
                    break;
                }
            }

            // ... if not, create one
            if (consensus_edge == nullptr) {
                consensus_edge = new_compact_edge(lcursor.cursor, 
                                                  rcursor.cursor,
                                                  edge_meta,
                                                  edge_sequence);
            } else {
                // if so, remove it from the set of found edges 
                // so we can resolve them
                segment_edges.erase(consensus_tag_edge_pair);
            }

            // remove all of the tags for this segment from stale edges
            // and map them to the consensus edge
            for (auto segment_tag: segment_tags) {

                for (auto tag_edge_pair: segment_edges) {
                    CompactEdge* stale_edge = tag_edge_pair.second;
                    stale_edge->tags.erase(segment_tag);
                }

                consensus_edge->tags.insert((HashIntoType)segment_tag);
                tags_to_edges[(HashIntoType)segment_tag] = consensus_edge;
                // don't forget to remove from the unresolved set
                unresolved_tags.erase(segment_tag);
            }

            // one last run through to delete stale edges and add remaining tags
            // to the unresolved set
            for (auto tag_edge_pair: segment_edges) {
                CompactEdge* edge = tag_edge_pair.second;
                if (edge != nullptr) {
                    for (auto stale_tag: edge->tags) {
                        unresolved_tags.insert(tag_edge_pair.first);
                        tags_to_edges.erase(stale_tag);
                    }
                }
                // don't forget to clean up
                delete_compact_edge(edge);
            }

            if (right_node != nullptr) {
                right_node->add_in_edge(
                        edge_sequence[edge_sequence.length()-(graph->ksize())],
                        consensus_edge);
            }

            if (left_node != nullptr) {
                left_node->add_out_edge(
                        edge_sequence[graph->ksize()],
                        consensus_edge);
            }

        }
    }

    void update_compact_dbg(const std::string& sequence) {
        std::cout << "update_compact_dbg()" << std::endl;

        UHashSet new_tags;
        UHashSet known_tags;
        uint64_t n_consumed = 0;
        graph->consume_sequence(sequence);
        consume_sequence_and_tag(sequence, n_consumed, new_tags, known_tags);
        update_compact_dbg(new_tags);
    }

    void index_sequence_tags(const std::string& seq,
                             KmerQueue unresolved_tags) {
        bool kmer_tagged;

        KmerIterator kmers(seq.c_str(), ksize());
        Traverser traverser(graph.get());
        Kmer kmer;

        KmerFilter finder = [&] (const Kmer& node) {
            
        };


        uint32_t since = tag_density / 2 + 1;

        while(!kmers.done()) {
            kmer = kmers.next();
            kmer_tagged = set_contains(tags_to_edges, kmer);
            if (kmer_tagged) {
                since = 1;
            } else {
                ++since;
            }

            if (since >= tag_density) {
                new_tags.insert(kmer);
                since = 1;
            }

            if (traverser.degree

        } // iteration over kmers

        if (since >= tag_density/2 - 1) {
            new_tags.insert(kmer);	// insert the last k-mer, too.
        }
    }

};



}




#endif
