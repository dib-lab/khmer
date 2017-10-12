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

class Junction {
protected:
    
    static uint64_t junction_counter;

public:

    // [u]-->[v (HDN)]-->[w]
    HashIntoType u;
    HashIntoType v;
    HashIntoType w;
    const uint64_t ID;

    Junction() = default;

    Junction(HashIntoType u, HashIntoType v, HashIntoType w) :
        u(u), v(v), w(w), ID(junction_counter++) {}

    uint64_t index() const { return u ^ v ^ w; }

    bool matches(HashIntoType u,
                 HashIntoType v,
                 HashintoType w) { return (u ^ v ^ w) == index(); }

    friend std::ostream& operator<< (std::ostream& stream,
                                    const Junction& j);
    friend bool operator== (const Junction& lhs,
                            const Junction& rhs) {
        return lhs.index() == rhs.index();
    }
};


typedef std::list<Junction*> JunctionList;

#define FW 1
#define RC 0

class LinkTreeSegment {
protected:
    static uint64_t segment_counter;
public:
    const uint64_t segment_id;
    std::vector<uint64_t> 
    std::vector<uint64_t> junctions; // IDs of junctions in this segment
    uint64_t children[4]; // child segment IDs
    uint32_t count;

    LinkTreeSegment() : segment_id(segment_counter++) {}

};



// Compare the two link cursors first by their age within the
// traversal, then by their global age. We prefer links
// that were found earlier in the traversal, but were created
// most recently.
/*
inline bool CompareLinks(const LinkCursor& a, const LinkCursor& b)
{
    if (a.traversal_age < b.traversal_age) { return false; }
    if (b.traversal_age < a.traversal_age) { return true; }

    if ((a.link)->time_created > (b.link)->time_created) { return false; }
    if ((b.link)->time_created > (a.link)->time_created) { return true; }

    return false;
}
*/


typedef std::unordered_multimap<HashIntoType, LinkSegment*> LinkSegmentMap;
typedef std::vector<LinkTreeSegment> LinkSegmentVector;
typedef std::vector<Junction> JunctionVector;
typedef std::unordered_map<HashIntoType, Junction*> JunctionMap;
typedef std::pair<HashIntoType, LinkSegment*> LinkMapPair;

class GraphLinker
{
protected:
    // map from starting high degree nodes to associated links
    LinkSegmentMap link_segments;
    // linear storage for Junctions
    JunctionVector junctions;
    // map from junction keys to Junction*
    JunctionMap    junction_map;
    uint64_t n_sequences_added;

    // can cause duplicate junctions if we don't check the map first;
    // this just keeps junctions on the same cache line if they're
    // allocated one after the other
    Junction* new_junction(HashIntoType u, HashIntoType v, HashIntoType w)
    {
        junctions.emplace_back(u, v, w);

        Junction* j = &(junctions.back);
        junctions[j->index()] = j;
        return j;
    }

    LinkTreeSegment* new_link_tree_segment(uint64_t junction_id)
    {
        link_segments.emplace_back();
        l = &(link_segments.back);
        return l;
    }

public:

    std::shared_ptr<Hashgraph> graph;
    
    GraphLinker(std::shared_ptr<Hashgraph> graph) :
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

    Junction* get_junction_by_index(HashIntoType key) const
    {
        auto search = junctions.find(key);
        if (search != junctions.end()) {
            return search->second;
        }
        return nullptr;
    }

    Junction* get_junction_by_index(HashIntoType u, HashIntoType v, HashIntoType w) const
    {
        return get_junction(u ^ v ^ w);
    }

    Junction* fetch_or_new_junction(HashIntoType u, HashIntoType v, 
                                    HashIntoType w, uint64_t& counter)
    {
        Junction* j = get_junction_by_index(u, v, w);
        if (j != nullptr) {
            j->count = j->count + 1;
        } else {
            counter++;
            j = new_junction(u, v, w);
            j->count = 1;
        }
        return j;
    }

    std::vector<uint64_t> get_junctions_ids(const std::string& sequence) const
    {
        KmerIterator kmers(sequence.c_str(), graph->ksize());
        std::vector<uint64_t> junctions;

        Kmer u = kmers.next();
        if (kmers.done()) {
            return junctions;
        }
        Kmer v = kmers.next();
        if (kmers.done()) {
            return junctions;
        }
        Kmer w = kmers.next();

        Junction* j;

        while(!kmers.done()) {
            j = get_junction_by_index(u, v, w);
            if (j != nullptr) {
                junctions->push_back(j->ID);
            }
            u = v;
            v = w;
            w = kmers.next();
        }

        return junctions;
    }

    LinkSegment* get_link_segment(uint64_t start_junction_id) {
        auto saerch = link_starts.find(start_junction_id);
        if (search != link_starts.end()) {
            return search->second;
        }
        return nullptr;
    }

    LinkSegment* fetch_or_new_link_segment(uint64_t start_junction_id) {
        LinkSegment* l = get_link_segment(start_junction_id);
        if (l != nullptr) {

        } else {
            l = new LinkTreeSegment();
        }
    }
 
    void build_link(const std::string& sequence,
                     Link* &fw_link, Link* &rc_link)
    {
        std::cout << "build_link()" << std::endl;
        KmerIterator kmers(sequence.c_str(), graph->ksize());
        Kmer u, v, w;
        uint64_t d = 0;

        u = kmers.next();
        ++d;
        if (kmers.done()) {
            fw_link = rc_link = nullptr;
            return;
        }
        v = kmers.next();
        ++d;
        if (kmers.done()) {
            fw_link = rc_link = nullptr;
            return;
        }
        w = kmers.next();
        ++d;
    
        uint64_t n_new_fw = 0;
        Junction* start = nullptr;

        // find starting HDN that will index the Link
        Traverser traverser(graph.get());
        while(!kmers.done() && start == nullptr) {
            if (traverser.degree(v) > 2) {
                std::cout << "  - build_links: found lead HDN" << u << std::endl;
                start = fetch_or_new_junction(u, v, w, n_new_fw);
            }

            u = v;
            v = w;
            w = kmers.next();
            ++d;
        }





        if (fw_link->size() < 1) {
            std::cout << "  - build_links: (fw_link) no (new) junctions found." << std::endl;
            delete fw_link;
            fw_link = nullptr;
        }

        if (rc_link->size() < 1) {
            std::cout << "  - build_links: (rc_link) no (new) junctions found." << std::endl;
            delete rc_link;
            rc_link = nullptr;
        } else {
            rc_link->flanking_distance = d - last_rc_pos;
        }
    }


    void add_links(const std::string& sequence)
    {
        std::cout << "add_links('" << sequence << "')" << std::endl;
        Link * fw_link = nullptr;
        Link * rc_link = nullptr;
        n_sequences_added++;
        build_links(sequence, fw_link, rc_link);

        if (fw_link != nullptr) {
            std::cout << "  - add_links: insert found fw_link" << std::endl;
            Junction* start = fw_link->start_junction();
            std::cout << "    * start junction: " << &start << std::endl;
            links.insert(LinkMapPair(start->v, fw_link));
        }
        if (rc_link != nullptr) {
            std::cout << "  - add_links: insert found rc_link" << std::endl;
            Junction* start = rc_link->start_junction();
            std::cout << "    * start junction: " << &start << std::endl;
            links.insert(LinkMapPair(start->v, rc_link));
        }
    }

    uint64_t get_links(Kmer hdn, std::shared_ptr<LinkList> found_links)
    {
        auto range = links.equal_range(hdn);
        uint64_t n_links_found = 0;
        for (auto it = range.first; it != range.second; ++it) {
            found_links->push_back(it->second);
            n_links_found++;
        }
        return n_links_found;
    }

    uint64_t get_links(std::list<Kmer> high_degree_nodes, 
                   std::shared_ptr<LinkList> found_links)
    {
        uint64_t n_links_found = 0;
        for (Kmer hdn : high_degree_nodes) {
            n_links_found += get_links(hdn, found_links);
        }
        return n_links_found;
    }

    std::shared_ptr<LinkList> get_links(const std::string& sequence)
    {
        std::cout << "get_links('" << sequence << "')" << std::endl;
        KmerIterator kmers(sequence.c_str(), graph->ksize());
        Traverser traverser(graph.get());
        std::list<Kmer> hdns;

        std::cout << "  - get_link: start iterating k-mers" << std::endl;
        Kmer kmer = kmers.next();
        while(!kmers.done()) {
            if (traverser.degree(kmer) > 2) {
                hdns.push_back(kmer);
            }
            kmer = kmers.next();
        }
        std::cout << "  - get_links: found " << hdns.size() << " hdns" << std::endl;

        std::shared_ptr<LinkList> found_links = make_shared<LinkList>();
        get_links(hdns, found_links);
        std::cout << "  - get_links: returning " << found_links->size() << " links" << std::endl;
        return found_links;
    }

};

/*

class LinkedAssembler
{
    std::shared_ptr<LinearAssembler> linear_asm;

public:

    const std::shared_ptr<Hashgraph> graph;
    const std::shared_ptr<GraphLinker> linker;
    WordLength _ksize;

    explicit LinkedAssembler(std::shared_ptr<GraphLinker> linker) :
        _ksize(linker->ksize()), graph(linker->graph), linker(linker)
    {
        linear_asm = make_shared<LinearAssembler>(graph.get());
    }

    StringVector assemble(const Kmer seed_kmer,
                          std::shared_ptr<Hashgraph> stop_bf=nullptr) const;

    template <bool direction>
    void _assemble_directed(AssemblerTraverser<direction>& cursor,
                            StringVector& paths) const;
};
*/


}




#endif
