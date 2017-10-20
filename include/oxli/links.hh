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
typedef std::unordered_set<HashIntoType> HashSet;
typedef std::vector<HashIntoType> HashVector;
typedef std::unordered_multimap<HashIntoType, uint64_t> HashIDMap;

class CompactEdge;
class CompactNode {
protected:
    static uint64_t node_counter;
public:
    Kmer kmer;
    uint32_t count;
    const uint64_t node_id;
    CompactEdge* in_edges[4] = {nullptr, nullptr, nullptr, nullptr};
    CompactEdge* out_edges[4] = {nullptr, nullptr, nullptr, nullptr};

    CompactNode(Kmer kmer): kmer(kmer), count(0), node_id(node_counter++) {}

    friend bool operator== (const CompactNode* lhs, const CompactNode* rhs) {
        return lhs->node_id == rhs->node_id;
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
};


class CompactEdge {
public:
    HashIntoType left; // left and right HDNs
    HashIntoType right;
    HashSet tags;
    bool is_tip;
    bool dirty;
    std::string sequence; // optional

    CompactEdge(HashIntoType left, HashIntoType right) : 
        left(left), right(right), is_tip(false), dirty(false) {}
    CompactEdge(HashIntoType left, HashIntoType right, bool is_tip) :
        left(left), right(right), is_tip(is_tip), dirty(false) {}
};


typedef std::vector<CompactNode> CompactNodeVector;
typedef std::unordered_multimap<HashIntpType, CompactEdge*> TagEdgeMap;
typedef std::pair<HashIntoType, CompactEdge*> TagEdgePair;
typedef std::set<TagEdgePair> TagEdgePairSet;


class GraphLinker
{
protected:
    // map from HDN hashes to CompactNode IDs
    HashIDMap hdn_ids;
    // linear storage for CompactNodes
    CompactNodeVector compact_nodes;
    // map from tags to CompactEdges
    TagEdgeMap tags_to_edges;

    uint64_t n_sequences_added;

    CompactNode* new_compact_node(Kmer hdn)
    {
        CompactNode * v = get_compact_node_by_kmer(hdn);
        if (v == nullptr) {
            compact_nodes.emplace_back(hdn);
            v = &(compact_nodes.back);
            hdn_ids[hdn] = v->node_id;
        }
        return v;
    }

public:

    shared_ptr<Hashgraph> graph;
    
    GraphLinker(shared_ptr<Hashgraph> graph) :
        graph(graph), n_sequences_added(0)
    {

    }

    WordLength ksize() const
    {
        return graph->ksize();
    }

    void report() const
    {
        std::cout << "GraphLinker (@" << this << " with "
            << "Hashgraph @" << graph.get() << ")" << std::endl;
        std::cout << "  * " << junctions.size() << " junctions" << std::endl;
        std::cout << "  * " << n_sequences_added << " sequences added" << std::endl;
        std::cout << "  * " << links.size() << " links" << std::endl;
    }

    CompactNode* get_compact_node_by_kmer(HashIntoType hdn) const
    {
        auto search = hdn_ids.find(hdn);
        if (search != hdn_ids.end()) {
            uint64_t ID = search->second;
            return &(compact_nodes[ID]);
        }
        return nullptr;
    }

    CompactNode* fetch_or_new_compact_node(Kmer hdn)
    {
        CompactNode* v = get_compact_node_by_kmer(hdn);
        if (v != nullptr) {
            v->count += 1;
        } else {
            v = new_compact_node(hdn);
            v->count = 1;
        }
        return v;
    }

    std::vector<uint64_t> get_compact_node_ids(const std::string& sequence) const
    {
        KmerIterator kmers(sequence.c_str(), graph->ksize());
        std::vector<uint64_t> ids;

        CompactNode* node;

        while(!kmers.done()) {
            Kmer kmer = kmers.next();

            node = get_compact_node_by_kmer(kmer);
            if (j != nullptr) {
                ids->push_back(node->node_id);
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

    TagEdgePair get_tag_edge_pair(uint64_t tag) {
        auto search = tags_to_edges.find(tag);
        if (search != tags_to_edges.end()) {
            return &search;
        }
        return nullptr;
    }

    CompactEdge* get_edge_froms_tags(HashSet& tags) {
        CompactEdge * edge = nullptr;
        for (auto tag: tags) {
            edge = get_compact_edge(tag);
            if (edge != nullptr) {
                break;
            }
        }
        return edge;
    }

    void assign_tags_to_edge(HashSet& tags, CompactEdge* edge) {
        for (auto tag: tags) {
            tags_to_edges[tag] = edge;
        }
    }

    KmerFilter get_tag_finder(shared_ptr<TagEdgePairSet> found_pairs,
                              shared_ptr<KmerSet> found_unresolved_tags,
                              shared_ptr<KmerSet> all_unresolved_tags) {
        KmerFilter finder = [=] (const Kmer& node) {
            TagEdgePair edge_pair = get_tag_edge_pair(node);
            if (edge_pair != nullptr) {
                found_pairs->insert(edge_pair);
            }
            if (set_contains(all_unresolved_tags, node) {
                found_unresolved_tags->insert(node);
            }
            return false; // never filter node, just gather while we're looking
        }       
    }

    template<bool direction>
    void update_compact_dbg(shared_ptr<KmerSet> unresolved_tags) {

        HashVector unresolved_q;
        for (auto tag: unresolved_tags) {
            unresolved_q.push_back(tag);
        }

        CompactingAssembler at(graph.get());
        shared_ptr<SeenSet> visited = make_shared<SeenSet>();
        while(unresolved_q.size() > 0) {
            auto start_tag = unresolved_q.back();
            unresolved_q.pop_back();

            if (!set_contains(&unresolved_tags, start_tag)) {
                continue;
            }
            
            shared_ptr<TagEdgePairSet> segment_edges = make_shared<TagEdgePairSet>();
            shared_ptr<HashSet> segment_new_tags = make_shared<HashSet>();
            CompactingAT<LEFT> lcursor(graph.get(), tag


        }

        shared_ptr<TagEdgePairSet> found_edges;
        shared_ptr<
        filters.push_back(get_tag_finder(found_edges));
        CompactingAT<direction> cursor(graph, from, visited);
        at._assemble_directed<direction>(cursor);

        if (cursor.cursor_degree() == 0) { // at a tip
            
        }
        
    }

    void update_compact_dbg(const std::string& sequence) {
        std::cout << "update_compact_dbg()" << std::endl;

        shared_ptr<HashSet> new_tags = make_shared<HashSet>();
        shared_ptr<HashSet> known_tags = make_shared<HashSet>();
        uint64_t n_consumed = 0;
        consume_sequence_and_tag(sequence, n_consumed, new_tags, known_tags);


    }

    // sigh at duplicated code from Hashgraph... perhaps in the future
    // that function can be templated to accept multiple Set types for tags
    void consume_sequence_and_tag(const std::string& seq,
                                  unsigned long long& n_consumed,
                                  shared_ptr<KmerSet> new_tags,
                                  shared_ptr<KmerSet> known_tags) {
        bool kmer_tagged;

        KmerIterator kmers(seq.c_str(), _ksize);
        Kmer kmer;

        unsigned int since = _tag_density / 2 + 1;

        while(!kmers.done()) {
            kmer = kmers.next();
            bool is_new_kmer;

            if ((is_new_kmer = graph->add(kmer))) {
                ++n_consumed;
            }

            if (is_new_kmer) {
                ++since;
            } else {
                kmer_tagged = set_contains(tags_to_edges, kmer);
                if (kmer_tagged) {
                    since = 1;
                    known_tags->insert(kmer);
                } else {
                    ++since;
                }
            }

            if (since >= _tag_density) {
                new_tags->insert(kmer);
                since = 1;
            }

        } // iteration over kmers

        if (since >= _tag_density/2 - 1) {
            new_tags->insert(kmer);	// insert the last k-mer, too.
        }
    }

};



}




#endif
