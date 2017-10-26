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
    Kmer out_hash;
    UHashSet tags;
    compact_edge_meta_t meta;
    std::string sequence;

    CompactEdge(Kmer in_node, Kmer out_node) : 
        in_node(in_node), out_node(out_node), meta(IS_FULL_EDGE) {}
    
    CompactEdge(Kmer in_node, Kmer out_node, compact_edge_meta_t meta) :
        CompactEdge(in_node, out_node), meta(meta), {}

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

    KmerFilter compact_node_filter;

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

        compact_node_filter = [&] (const Kmer& node) {
            return get_compact_node_by_kmer(node) != nullptr;
        };
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
    }

    KmerHelper get_tag_collector(KmerSet& tags) {
        KmerHelper collector = [&] (const Kmer& node) {
            if (set_contains(tags_to_egdes, node)) {
                tags.insert(node);
            }
        };
    }

    void update_compact_dbg(KmerQueue& unresolved_q) {

        KmerIterator kmers(sequence.c_str(), ksize());

        KmerSet induced_hdns;
        Kmer kmer = kmers.next();
        kmer.set_forward();
        CompactingAT<TRAVERSAL_LEFT> lcursor(graph.get(), kmer);
        CompactingAT<TRAVERSAL_RIGHT> rcursor(graph.get(), kmer);

        while(!kmers.done()) {
            if(lcursor.degree(kmer) > 1 || rcursor.degree(kmer) > 1) {
                if (fetch_or_new_compact_node(kmer)->count == 1) {
                    induced_hdns.insert(kmer);
                }
            }
            kmer = kmers.next();
        }
        CompactNode* left_flank=nullptr, right_flank=nullptr;
        if ((left_flank = get_compact_node_by_kmer(lcursor.cursor)) != nullptr) {
            
        }


        while(!induced_hdns.empty()) {
            Kmer root_kmer = *induced_hdns.begin();
            induced_hdns.erase(root_hdn);
            root_kmer.set_forward();

            CompactNode* root_node = get_compact_node_by_kmer(root_kmer);

            KmerQueue left_neighbors;
            lcursor.neighbors(root_kmer, left_neighbors);
            while(!left_neighbors.empty()) {
                Kmer neighbor = left_neighbors.back();
                left_neighbors.pop_back();
                lcursor.set_cursor(neighbor);

                KmerSet tags_left;
                lcursor.push_helper(get_tag_collector(tags_left));
                std::string segment_seq = at._assemble_directed(lcursor);
                
                CompactNode* left_node = get_compact_node_by_kmer(lcursor.cursor);
                CompactEdge* segment_edge = root_node->get_in_edge(*segment_seq.rbegin());

                if (segment_edge == nullptr) {
                    if (!tags_left.empty()) {
                        segment_edge = get_compact_edge(*tags_left.begin());
                        
                    }
                }


                } else {
                    segment_edge =
                }
                    if (segment_edge->in_node == lcursor.cursor) {
                        // so far so good
                        if (segment_edge->out_node == root_node->kmer) {
                            // this edge is fine

                        } else {
                            // root must have split this edge
                            // this edge's old out node must be right of root
                            // let's fix it!
                            
                            for (auto edge_tag : segment_edge->tags) {
                                if (!set_contains(segment_tags, edge_tag)) {
                                    tags_to_edges.erase(edge_tag);
                                }
                            }
                        }
                    } else {
                        
                    }
                }

            }

            lcursor.set_cursor(root_hdn);
            rcursor.set_cursor(root_hdn);

            KmerSet tags_left, tags_right;

            rcursor.push_helper(get_tag_collector(tags_right));

            std::string left_seq = at._assemble_directed(lcursor);
            std::string right_seq = at._assemble_directed(rcursor);



        }



        CompactingAssembler at(graph.get());
        CompactingAT<TRAVERSAL_LEFT> lcursor(graph.get(), seed_kmer,
                                                 filters);
        CompactingAT<TRAVERSAL_RIGHT> rcursor(graph.get(), seed_kmer,
                                                  filters);
        while(unresolved_q.size() > 0) {
            HashIntoType seed_tag = unresolved_q.back();
            Kmer seed_kmer = graph->build_kmer(seed_tag);
            seed_kmer.set_forward();
            unresolved_q.pop_back();
            pdebug("iter unresolved q, seed=" << seed_kmer << 
                   ", " << unresolved_tags.size() << " unresolved tags left");
            if (!set_contains(unresolved_tags, seed_tag)) {
                pdebug("skip tag not in unresolved set");
                continue;
            }

            // build filters to prepare for traversal...
            // first, tag gatherer
            CompactEdgeSet segment_edges;
            UHashSet segment_tags;
            KmerFilterList filters;
            filters.push_back(get_tag_finder(segment_edges,
                                             segment_tags,
                                             unresolved_tags));

            // ... and a shared list of known k-mers
            shared_ptr<SeenSet> visited = make_shared<SeenSet>();
            filters.push_back(get_visited_filter(visited));
            filters.push_back(compact_node_filter);

            CompactingAT<TRAVERSAL_LEFT> lcursor(graph.get(), seed_kmer,
                                                 filters);
            CompactingAT<TRAVERSAL_RIGHT> rcursor(graph.get(), seed_kmer,
                                                  filters);

            // assemble both dirs until HDN, collecting tags along the way
            std::string left_seq = at._assemble_directed(lcursor);
            std::string right_seq = at._assemble_directed(rcursor);
            std::string edge_sequence = left_seq + right_seq.substr(graph->ksize());

            pdebug("assembled linear path: " << segment_edges.size() <<
                   " existing edges, " << segment_tags.size() << 
                   " tags gathered, sequence=" << edge_sequence);

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
                pdebug("has right HDN=" << *right_node);

            }
            // same for other side
            KmerQueue lneighbors;
            if (lcursor.neighbors(lneighbors) > 1 || 
                rcursor.neighbors(lcursor.cursor, lneighbors) > 1) {

                left_node = fetch_or_new_compact_node(lcursor.cursor);
                pdebug("has left HDN=" << *left_node);
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
            for (auto intersecting_edge: segment_edges) {
                if (intersecting_edge->meta != edge_meta) {
                    continue;
                }

                if (intersecting_edge->in_hash == (HashIntoType)lcursor.cursor &&
                    intersecting_edge->out_hash == (HashIntoType)rcursor.cursor) {

                    consensus_edge = intersecting_edge;
                    break;
                }
            }

            // ... if not, create one
            if (consensus_edge == nullptr) {
                pdebug("create consensus edge");
                consensus_edge = new_compact_edge(lcursor.cursor, 
                                                  rcursor.cursor,
                                                  edge_meta,
                                                  edge_sequence);
            } else {
                // if so, remove it from the set of found edges 
                // so we can resolve them
                pdebug("use existing consensus edge");
                segment_edges.erase(consensus_edge);
            }

            // remove all of the tags for this segment from stale edges
            // and map them to the consensus edge
            for (auto segment_tag: segment_tags) {

                for (auto stale_edge: segment_edges) {
                    stale_edge->tags.erase(segment_tag);
                }

                consensus_edge->tags.insert((HashIntoType)segment_tag);
                tags_to_edges[(HashIntoType)segment_tag] = consensus_edge;
                // don't forget to remove from the unresolved set
                unresolved_tags.erase(segment_tag);
            }

            // one last run through to delete stale edges and add remaining tags
            // to the unresolved set
            for (auto stale_edge: segment_edges) {
                if (stale_edge != nullptr) {
                    for (auto stale_tag: stale_edge->tags) {
                        unresolved_tags.insert(stale_tag);
                        tags_to_edges.erase(stale_tag);
                    }
                }
                // don't forget to clean up
                delete_compact_edge(stale_edge);
            }

            if (right_node != nullptr) {
                right_node->add_in_edge(
                        edge_sequence[edge_sequence.length()-(graph->ksize())],
                        consensus_edge);
                if (right_node->count == 1) {
                    while(rneighbors.size() > 0) {
                        // if it was a new HDN, we it could be between tag and
                        // nearest existing HDN; so, we have to check from its
                        // neighbors
                        auto rn = rneighbors.back();
                        rneighbors.pop_back();
                        unresolved_q.push_back(rn);
                        unresolved_tags.insert(rn);
                    }
                }
            }

            if (left_node != nullptr) {
                left_node->add_out_edge(
                        edge_sequence[graph->ksize()],
                        consensus_edge);
                if (left_node->count == 1) {
                    while(lneighbors.size() > 0) {
                        auto ln = lneighbors.back();
                        lneighbors.pop_back();
                        unresolved_q.push_back(ln);
                        unresolved_tags.insert(ln);
                    }
                }
            }

        }
    }

    void update_compact_dbg(const std::string& sequence) {
        pdebug("sequence=" << sequence);

        KmerQueue unresolved_tags;
        uint64_t prev_n_kmers = graph->n_unique_kmers();
        graph->consume_string(sequence);
        if (graph->n_unique_kmers() - prev_n_kmers) {
            seed_sequence_tags(sequence, unresolved_tags);
            update_compact_dbg(unresolved_tags);
        }
    }

    void orient_sequence(const std::string& seq,
                             KmerQueue& unresolved_tags) {

        KmerIterator kmers(seq.c_str(), ksize());
        Traverser traverser(graph.get());
        Kmer kmer = kmer.next();

        

        while(!kmers.done()) {
            kmer = kmers.next();
            if (set_contains(tags_to_edges, kmer)) {
                since = 1;
            } else {
                ++since;
            }

            if (since >= tag_density) {

                if(!compact_node_filter(kmer)) {
                    unresolved_tags.push_back(kmer);
                }

                since = 1;
            }

            if (traverser.degree_left(kmer) > 1 ||
                traverser.degree_right(kmer) > 1) {
                
                pdebug("found HDN " << kmer);
                CompactNode* hdn = get_compact_node_by_kmer(kmer);
                if (hdn == nullptr) { // new HDN
                    pdebug("HDN is new, seed neighbors. " << 
                            unresolved_tags.size() << " current tags");

                    traverser.traverse(kmer, unresolved_tags);

                    pdebug(unresolved_tags.size() << " after seeding");
                }

            }

        } // iteration over kmers

        if (since >= tag_density/2 - 1) {
            if (!compact_node_filter(kmer)) {
                unresolved_tags.push_back(kmer);	// insert the last k-mer, too.
            }
        }

        pdebug("seeded " << unresolved_tags.size() << " tags");
    }

};



}


#endif
