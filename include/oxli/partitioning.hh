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

#include "gmap.hh"
#include "oxli.hh"
#include "kmer_hash.hh"
#include "hashtable.hh"
#include "hashgraph.hh"
#include "kmer_filters.hh"
#include "traversal.hh"

#ifndef DEBUG_SP
#define DEBUG_SP 0
#endif

namespace oxli
{


class Component;
typedef std::shared_ptr<Component> ComponentPtr;


class ComponentPtrCompare {
    public:
        bool operator() (const ComponentPtr& lhs, const ComponentPtr& rhs) const;
};


class ComponentPtrCompare;
typedef std::set<ComponentPtr, ComponentPtrCompare> ComponentPtrSet;
typedef std::vector<ComponentPtr> ComponentPtrVector;
typedef std::vector<HashIntoType> TagVector;
typedef GuardedHashMap<ComponentPtr, false> GuardedHashCompMap;


class Component {

    private:
        
        static uint64_t n_created;
        static uint64_t n_destroyed;
        bool alive;

    public:

        const uint64_t component_id;
        TagVector tags;

        explicit Component(): component_id(n_created), alive(true) {
            n_created++;
        }

        explicit Component(uint64_t component_id): component_id(component_id) {
            n_created++;
        }

        ~Component() {
            n_destroyed++;
        }

        void kill() {
            tags.clear();
            alive = false;
        }

        bool is_alive() const {
            return alive;
        }

        uint64_t get_n_created() const {
            return n_created;
        }

        uint64_t get_n_destroyed() const {
            return n_destroyed;
        }

        void add_tag(HashIntoType tag) {
            tags.push_back(tag);
        }

        void add_tags(TagVector& new_tags) {
            tags.insert(tags.end(),
                        new_tags.begin(),
                        new_tags.end());
        }

        uint64_t get_n_tags() const {
            return tags.size();
        }

        friend bool operator==(const Component& lhs, const Component& rhs) {
            return lhs.component_id == rhs.component_id;
        }

        friend bool operator<(const Component& lhs, const Component& rhs) {
            return lhs.component_id < rhs.component_id;
        }

        friend std::ostream& operator<< (std::ostream& stream, const Component& comp);
};


class ComponentMap {

    private:
    
        // We should exclusively own tag_component_map.
        std::shared_ptr<GuardedHashCompMap> tag_component_map;
        std::shared_ptr<ComponentPtrVector> components;
        uint32_t components_lock;
        uint64_t component_counter;
        uint64_t n_live_components;

    public:


        explicit ComponentMap(WordLength ksize,
                              WordLength n_tables,
                              uint64_t max_table_size);

        void create_component(TagVector& tags);
        uint32_t create_and_merge_components(TagVector& tags);
        void map_tags_to_component(TagVector& tags, ComponentPtr& comp);
        uint32_t merge_components(ComponentPtr& root, ComponentPtrSet& comps);

        bool contains(HashIntoType tag)
        {
            return tag_component_map->contains(tag);
        }

        ComponentPtr get(HashIntoType tag) const {
            return tag_component_map->get(tag);
        }

        uint64_t get_n_components() const {
            return n_live_components;
        }

        uint64_t get_n_tags() const {
            return tag_component_map->size();
        }

        std::weak_ptr<ComponentPtrVector> get_components() const {
            return std::weak_ptr<ComponentPtrVector>(components);
        }

        std::weak_ptr<GuardedHashCompMap> get_tag_component_map() const {
            return std::weak_ptr<GuardedHashCompMap>(tag_component_map);
        }

        inline void acquire_components() {
            while(!__sync_bool_compare_and_swap( &components_lock, 0, 1));
        }

        inline void release_components() {
            __sync_bool_compare_and_swap( &components_lock, 1, 0);
        }
};


class StreamingPartitioner {

    private:
    
        uint32_t _tag_density;
        std::shared_ptr<ComponentMap> partitions;
        uint64_t n_consumed;

    public:
        // We're not graph's owner, simply an observer.
        // Unforunately our ownership policies elsewhere are a mess
        Hashgraph * graph;
        //std::weak_ptr<Hashgraph> graph;

        explicit StreamingPartitioner(Hashgraph * graph, 
                                      uint32_t tag_density=DEFAULT_TAG_DENSITY);

        uint64_t consume(const std::string& seq);
        uint64_t consume_pair(const std::string& first,
                          const std::string& second);

        uint64_t seed_sequence(const std::string& seq,
                              TagVector& tags,
                              KmerQueue& seeds,
                              std::set<HashIntoType>& seen);

        uint64_t consume_fasta(std::string const &filename);
        void find_connected_tags(KmerQueue& node_q,
                                 TagVector& found_tags,
                                 std::set<HashIntoType>& seen,
                                 bool truncate=false) const;

        ComponentPtr get_tag_component(std::string& kmer) const;
        ComponentPtr get_nearest_component(Kmer kmer) const;
        ComponentPtr get_nearest_component(std::string& kmer) const;

        std::weak_ptr<ComponentPtrVector> get_components() const {
            return partitions->get_components();
        }

        std::weak_ptr<GuardedHashCompMap> get_tag_component_map() const {
            return partitions->get_tag_component_map();
        }

        uint64_t get_n_components() const {
            return partitions->get_n_components();
        }

        uint64_t get_n_tags() const {
            return partitions->get_n_tags();
        }

        uint32_t get_tag_density() const {
            return _tag_density;
        }
};


}

#endif
