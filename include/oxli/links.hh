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


namespace oxli {

struct Junction {
    HashIntoType u;
    HashIntoType v;
    //uint32_t distance_prev;
    uint64_t count;
    Junction() = default;
    HashIntoType id() const { return u ^ v; }


    friend std::ostream& operator<< (std::ostream& stream,
                                    const Junction& j);
};


typedef std::list<Junction*> JunctionList;


class Link {
private:

    static uint64_t n_links;

protected:

    JunctionList junctions;
    bool forward;
    uint64_t link_id;
    uint64_t time_created; // in "read space" ie read number
    //size_t extent;

public:
    
    Link(uint64_t time_created, bool forward=true) :
        forward(forward), link_id(n_links),
        time_created(time_created)
    {
        n_links++;
    }
    inline void push_back(Junction* junction)
    {
        //if (junctions.size() == 0) {
        //    junction->start = true;
            //junction->distance_prev = 0;
        //}
        junctions.push_back(junction);
    }

    inline void push_front(Junction* junction)
    {
        //if (junctions.size() > 0) {
        //    start_junction()->start = false;
        //}
        //junction->start = true; // just to be sure
        //junction->distance_prev = 0;
        junctions.push_front(junction);
    }

    void insert_junction(JunctionList::iterator insert_before,
                         Junction* junction)
    {
        if (insert_before == begin()) {
            push_front(junction);
        } else {
            junctions.insert(insert_before, junction);
        }
    }

    inline Junction* start_junction() const
    {
        return junctions.front();
    }

    inline Junction* end_junction() const
    {
        return junctions.back();
    }

    inline bool is_forward() const
    {   
        return forward;
    }

    inline size_t size() const
    {
        return junctions.size();
    }

    inline JunctionList::iterator begin() 
    {
        return junctions.begin();
    }

    inline JunctionList::iterator end() 
    {
        return junctions.end();
    }

    JunctionList& get_junctions() {
        return junctions;
    }

};

typedef std::unordered_multimap<HashIntoType, Link*> LinkMap;
typedef std::list<Link*> LinkList;
typedef std::unordered_map<HashIntoType, Junction*> JunctionMap;
typedef std::pair<HashIntoType, Link*> LinkMapPair;

class GraphLinker
{
protected:
    // map from starting high degree nodes to associated links
    LinkMap links;
    // map from junction keys to Junction*
    JunctionMap junctions;
    std::shared_ptr<Hashgraph> graph;
    uint64_t n_sequences_added;

public:
    
    GraphLinker(std::shared_ptr<Hashgraph> graph) :
        graph(graph), n_sequences_added(0)
    {

    }

    void report() const
    {
        std::cout << "GraphLinker (@" << this << " with "
            << "Hashgraph @" << graph.get() << ")" << std::endl;
        std::cout << "  * " << junctions.size() << " junctions" << std::endl;
        std::cout << "  * " << n_sequences_added << " sequences added" << std::endl;
        std::cout << "  * " << links.size() << " links" << std::endl;
    }

    Junction* get_junction(HashIntoType key) const
    {
        auto search = junctions.find(key);
        if (search != junctions.end()) {
            return search->second;
        }
        return nullptr;
    }

    Junction* get_junction(HashIntoType u, HashIntoType v) const
    {
        return get_junction(u ^ v);
    }

    Junction* get_junction(Junction& junction) const
    {
        return get_junction(junction.id());
    }

    Junction* new_junction(HashIntoType u, HashIntoType v)
    {
        Junction* j = new Junction();
        j->u = u;
        j->v = v;
        j->count = 0;
        junctions[j->id()] = j;
        return j;
    }

    Junction* fetch_or_new_junction(HashIntoType u, HashIntoType v, uint64_t& counter)
    {
        Junction* j = get_junction(u, v);
        if (j != nullptr) {
            j->count = j->count + 1;
        } else {
            counter++;
            j = new_junction(u, v);
            j->count = 1;
        }
        return j;
    }

    std::shared_ptr<JunctionList> get_junctions(const std::string& sequence) const
    {
        KmerIterator kmers(sequence.c_str(), graph->ksize());
        std::shared_ptr<JunctionList> junctions = make_shared<JunctionList>();

        Kmer u = kmers.next();
        if (kmers.done()) {
            return junctions;
        }
        Kmer v = kmers.next();

        Junction* j;

        while(!kmers.done()) {
            j = get_junction(u, v);
            if (j != nullptr) {
                junctions->push_back(j);
            }
            u = v;
            v = kmers.next();
        }

        return junctions;
    }
 
    void build_links(const std::string& sequence,
                     Link* &fw_link, Link* &rc_link)
    {
        std::cout << "build_links()" << std::endl;
        KmerIterator kmers(sequence.c_str(), graph->ksize());
        Kmer u, v;

        u = kmers.next();
        if (kmers.done()) {
            fw_link = rc_link = nullptr;
            return;
        }
        v = kmers.next();
    
        std::cout << "  - build_links: allocate new Link*" << std::endl;
        fw_link = new Link(n_sequences_added);
        rc_link = new Link(n_sequences_added, false);
        uint64_t n_new_fw = 0;
        uint64_t n_new_rc = 0;

        Traverser traverser(graph.get());
        while(!kmers.done()) {
            if (traverser.degree_right(u) > 1) {
                std::cout << "  - build_links: found FW HDN " << u << std::endl;
                fw_link->push_back(fetch_or_new_junction(u, v, n_new_fw));
            }
            if (traverser.degree_left(v) > 1) {
                std::cout << "  - build_links: found RC HDN " << v << std::endl;
                rc_link->push_front(fetch_or_new_junction(v, u, n_new_rc));
            }

            u = v;
            v = kmers.next();
        }

        if (fw_link->size() < 1 || n_new_fw == 0) {
            std::cout << "  - build_links: (fw_link) no (new) junctions found." << std::endl;
            delete fw_link;
            fw_link = nullptr;
        }

        if (rc_link->size() < 1 || n_new_rc == 0) {
            std::cout << "  - build_links: (rc_link) no (new) junctions found." << std::endl;
            delete rc_link;
            rc_link = nullptr;
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
            links.insert(LinkMapPair(start->u, fw_link));
        }
        if (rc_link != nullptr) {
            std::cout << "  - add_links: insert found rc_link" << std::endl;
            Junction* start = rc_link->start_junction();
            std::cout << "    * start junction: " << &start << std::endl;
            links.insert(LinkMapPair(start->u, rc_link));
        }
    }

    void get_links(std::list<HashIntoType> high_degree_nodes, 
                   std::shared_ptr<LinkList> found_links)
    {
        for (HashIntoType hdn : high_degree_nodes) {
            auto range = links.equal_range(hdn);
            for (auto it = range.first; it != range.second; ++it) {
                found_links->push_back(it->second);
            }
        }
    }

    std::shared_ptr<LinkList> get_links(const std::string& sequence)
    {
        std::cout << "get_links('" << sequence << "')" << std::endl;
        KmerIterator kmers(sequence.c_str(), graph->ksize());
        Traverser traverser(graph.get());
        std::list<HashIntoType> hdns;


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



}




#endif
