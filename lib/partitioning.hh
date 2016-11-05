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
#ifndef PARTITIONING_HH 
#define PARTITIONING_HH

#include <functional>
#include <memory>

#include "khmer.hh"
#include "kmer_hash.hh"
#include "hashtable.hh"
#include "hashbits.hh"
#include "counting.hh"
#include "kmer_filters.hh"
#include "traversal.hh"


namespace khmer
{

template <typename T>
class GuardedKmerMap {

    public:

        // Filter should be owned exclusively by GuardedKmerMap
        std::unique_ptr<Hashbits> filter;
        std::map<HashIntoType, T> data;
        
        explicit GuardedKmerMap(WordLength ksize,
                                unsigned short n_tables,
                                uint64_t max_table_size)
        {
            std::vector<uint64_t> table_sizes = get_n_primes_near_x(n_tables, max_table_size);
            filter = std::unique_ptr<Hashbits>(new Hashbits(ksize, table_sizes));
        }

        T get(HashIntoType kmer) {
            if (filter->get_count(kmer)) {
                auto search = data.find(kmer);
                if (search != data.end()) {
                    return search->second;
                }
            }
            
            return NULL;
        }

        void set(HashIntoType kmer, T item) {
            filter->count(kmer);
            data[kmer] = item;
        }

        bool contains(HashIntoType kmer) {
            return get(kmer) != NULL;
        }

};


class Component {

    private:
        
        static uint64_t n_created;
        uint64_t n_merges;

    public:

        const uint64_t component_id;
        std::set<HashIntoType> tags;

        explicit Component(): component_id(n_created), n_merges(0) {
            n_created++;
        }

        ~Component() {
            tags.clear(); // maybe not necessary?
        }

        void merge(std::set<std::shared_ptr<Component>> other_comps) {
            for (auto other : other_comps) {
                if (other.get() == this) {
                    continue;
                }
                this->add_tag(other->tags);
                this->n_merges = other->get_n_merges() + 1;
            }
        }

        void add_tag(HashIntoType tag) {
            tags.insert(tag);
        }

        void add_tag(std::set<HashIntoType> new_tags) {
            for (auto tag: new_tags) {
                add_tag(tag);
            }
        }

        uint64_t get_n_tags() {
            return tags.size();
        }

        uint64_t get_n_merges() {
            return n_merges;
        }
};
static uint64_t Component::n_created = 0;
typedef std::shared_ptr<Component> ComponentPtr;
typedef GuardedKmerMap<ComponentPtr> GuardedKmerCompMap;

class StreamingPartitioner {

    private:
    
        uint32_t _tag_density;

        // We're not graph's owner, simply an observer.
        std::weak_ptr<Hashtable> graph;
        // We should exlusively own tag_component_map.
        std::unique_ptr<GuardedKmerCompMap> tag_component_map;
        uint64_t n_components;

    public:

        explicit StreamingPartitioner(std::weak_ptr<Hashtable>& graph);

        void consume_sequence(const std::string& seq);
        void map_tags_to_component(std::set<HashIntoType> tags, ComponentPtr& comp);
        void find_connected_tags(KmerQueue& node_q,
                                 std::set<HashIntoType>& found_tags,
                                 std::set<HashIntoType>& seen);
        uint64_t get_n_components() {
            return n_components;
        }
};


}

#endif
