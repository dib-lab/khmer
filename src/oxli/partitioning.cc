#include <queue>
#include <set>
#include <deque>
#include <memory>
#include <functional>
#include <algorithm>
#include <iostream>

#include "oxli/hashtable.hh"
#include "oxli/hashgraph.hh"
#include "oxli/partitioning.hh"

using namespace oxli;
using namespace oxli::read_parsers;

uint64_t Component::n_created = 0;
uint64_t Component::n_destroyed = 0;

bool ComponentPtrCompare::operator() (const ComponentPtr& lhs, 
                                      const ComponentPtr& rhs) const {
    return *lhs < *rhs;
}

inline std::ostream& operator<< (std::ostream& stream, Component& comp) {
    stream << "<Component (id=" << comp.component_id << ", n_tags=" 
           << comp.get_n_tags() << ")>";
    return stream;
}

////////////////////////////////////////////////////////////////////////////////


ComponentMap::ComponentMap(WordLength ksize,
                           WordLength n_tables,
                           uint64_t max_table_size) : components_lock(0),
                                                      component_counter(0),
                                                      n_live_components(0)
{

    tag_component_map = std::unique_ptr<GuardedHashCompMap>(
                            new GuardedHashCompMap(ksize, 
                                                   n_tables, 
                                                   max_table_size));
    components = std::make_shared<ComponentPtrVector>();
}

void ComponentMap::map_tags_to_component(TagVector& tags,
                                         ComponentPtr& comp)
{
    for (auto tag: tags) {
        tag_component_map->set(tag, comp);
        comp->add_tag(tag);
    }
}

void ComponentMap::create_component(TagVector& tags)
{
    ComponentPtr new_comp = std::make_shared<Component>(component_counter);
    component_counter++;
    n_live_components++;
    components->push_back(new_comp);
    map_tags_to_component(tags, new_comp);

    //std::cout << "new component=" << *new_comp << std::endl;
    //std::cout << components->size() << " components in vector" << std::endl;
}


uint32_t ComponentMap::create_and_merge_components(TagVector& tags)
{

        // Now resolve components. First, get components from existing tags.
        ComponentPtrSet found_comps;
        TagVector new_tags;
        for (auto tag: tags) {
            ComponentPtr comp;
            if ((comp = tag_component_map->get(tag)) != NULL) {
                found_comps.insert(comp);
            } else {
                new_tags.push_back(tag);
            }
        }
        
        uint32_t n_merged = 1;
        if (found_comps.size() == 0) {
            create_component(tags);
        } else {
            // Choose the largest component as the root
            // We want to minimize tag copying
            ComponentPtr root_comp = *(found_comps.begin());
            for (auto other : found_comps) {
                if (other->get_n_tags() > root_comp->get_n_tags()) {
                    root_comp = other;
                }
            }
            // map the new tags to this component
            root_comp->add_tags(new_tags);
            map_tags_to_component(new_tags, root_comp);
            if (found_comps.size() > 1) {
                n_merged = merge_components(root_comp, found_comps);
            }
        }
        return n_merged;
}


uint32_t ComponentMap::merge_components(ComponentPtr& root, 
                                        ComponentPtrSet& comps)
{
    uint32_t n_merged = 1;
    //std::cout << "Merge with root=" << *root << std::endl;
    for (auto other : comps) {
        //std::cout << "\tmerge in " << *other << std::endl;
        if (*other == *root) {
            continue;
        }
        root->add_tags(other->tags); // transfer the tags from the other comp
        map_tags_to_component(other->tags, root);
        (*components)[other->component_id]->kill();
        (*components)[other->component_id] = nullptr;
        n_live_components--;
        n_merged++;

    }
                   // and active Python wrapper; this leaves them as sole owners
    return n_merged;
}


////////////////////////////////////////////////////////////////////////////////


StreamingPartitioner::StreamingPartitioner(Hashgraph * graph,
                                           uint32_t tag_density) :
    ComponentMap::ComponentMap(graph->ksize(),
                               graph->n_tables(),
                               _cstr_get_max_table_size(graph)),
    graph(graph),
    _tag_density(tag_density), 
    n_consumed(0)
{
}


uint64_t StreamingPartitioner::_cstr_get_max_table_size(Hashgraph * graph)
{
    std::vector<uint64_t> graph_table_sizes = graph->get_tablesizes(); 
    return  *std::max_element(graph_table_sizes.begin(),
                              graph_table_sizes.end());
}


uint64_t StreamingPartitioner::consume_fasta(const std::string& filename)
{
    ReadParserPtr<FastxReader> parser = get_parser<FastxReader>(filename);
    Read read;
    uint64_t n_consumed = 0;

    while (!parser->is_complete()) {
        if (n_consumed && (n_consumed % 10000 == 0)) {
            std::cout << "consumed " << n_consumed << "..." << std::endl;
        }
        try {
            read = parser->get_next_read( );
        } catch (NoMoreReadsAvailable) {
            break;
        }

        read.set_clean_seq();
        consume(read.sequence);
        n_consumed++;
    }

    return n_consumed;
}


uint64_t StreamingPartitioner::consume(const std::string& seq)
{
    TagVector tags;
    KmerQueue seeds;
    std::set<HashIntoType> seen;

    uint64_t n_new = seed_sequence(seq, tags, seeds, seen);
    find_connected_tags(seeds, tags, seen, false);
    //acquire_components();
    create_and_merge_components(tags);
    //release_components();
    return n_new;
}


uint64_t StreamingPartitioner::consume_pair(const std::string& first,
                                            const std::string& second)
{
    TagVector tags;
    KmerQueue seeds;
    std::set<HashIntoType> seen;

    uint64_t n_new = seed_sequence(first, tags, seeds, seen);
    n_new += seed_sequence(second, tags, seeds, seen);
    find_connected_tags(seeds, tags, seen, false);
    //acquire_components();
    create_and_merge_components(tags);
    //release_components();
    return n_new;
}


ComponentPtr StreamingPartitioner::get(std::string& kmer) const
{
    HashIntoType h = graph->hash_dna(kmer.c_str());
    return ComponentMap::get(h);
}


ComponentPtr StreamingPartitioner::get(HashIntoType h) const
{
    return ComponentMap::get(h);
}


uint64_t StreamingPartitioner::seed_sequence(const std::string& seq,
                                             TagVector& tags,
                                             KmerQueue& seeds,
                                             std::set<HashIntoType>& seen)
{
    /* For the following comments, let G be the set of k-mers
     * known in the graph before inserting the k-mers R from
     * &seq, with / as the difference, + as the union, and &
     * as the intersect operator.
     */
    //if (auto graphptr = graph.lock()) {
#if(SP_DEBUG)
    std::cout << "Consume sequence." << std::endl;
#endif
    uint64_t n_new = 0;
    ++n_consumed;

    if(graph != NULL) {
        KmerIterator kmers(seq.c_str(), graph->ksize());
        unsigned int since = _tag_density / 2 + 1;

        KmerSet intersection;

        bool in_known_territory = false;
        bool found_tag_in_territory = false;
        
        Kmer kmer;
        do {
            kmer = kmers.next();
            bool is_new_kmer = graph->add(kmer);
            bool kmer_tagged = false;

            if (is_new_kmer) {
                // A k-mer from U/G must be searched from for tags,
                // as it could be adjacent to a a k-mer in G/U
                if (in_known_territory && found_tag_in_territory) {
                    // If we had found a tag in the U&G component we just
                    // left, add the component to the seen set.
                    seen.insert(intersection.begin(), intersection.end());
                } /*else {
                    for (auto km : intersection) {
                        seeds.push(km);
                    }
                }*/
                intersection.clear();

                seeds.push_back(kmer);
                in_known_territory = false;
                found_tag_in_territory = false;
                ++since;
                ++n_new;
            } else {
                // Keep track of connected components in U&G: when we exit
                // this component, if there is a tag, we will want to add its nodes
                // to the seen set, as we do not need to traverse from them in the tag search.
                intersection.insert(kmer);
                in_known_territory = true;
                kmer_tagged = this->contains(kmer);
                if (kmer_tagged) {
                    since = 1;
                    tags.push_back(kmer);
                    found_tag_in_territory = true;
                } else {
                    ++since;
                }
            }

            if (since >= _tag_density) {
                tags.push_back(kmer);
                since = 1;
            }
        } while (!kmers.done());

        // always tag the last k-mer
        if (since >= _tag_density / 2) {
            tags.push_back(kmer);
        }
        seeds.push_back(kmer);

        // now go back and make sure to search from the first k-mer
        kmer = kmers.first();
        seeds.push_back(kmer);

#if(DEBUG_SP)
        std::cout << "Done iterating k-mers" << std::endl;
        std::cout << tags.size() << " tags in sequence" << std::endl;
#endif
    } else {
        throw oxli_ptr_exception("Hashgraph has been deleted.");
    }

    return n_new;
}

ComponentPtr StreamingPartitioner::find_nearest_component(std::string& kmer) const
{
    Kmer hashed = graph->build_kmer(kmer);
    return find_nearest_component(hashed);
}


ComponentPtr StreamingPartitioner::find_nearest_component(Kmer kmer) const
{
    TagVector tags;
    std::set<HashIntoType> seen;
    KmerQueue node_q;
    node_q.push_front(kmer);

    find_connected_tags(node_q, tags, seen, true);
    if (tags.size() > 0) {
        HashIntoType tag = *(tags.begin());
        return this->get(tag);
    } else {
        return NULL;
    }
}


void StreamingPartitioner::find_connected_tags(KmerQueue& node_q,
                                               TagVector& found_tags,
                                               std::set<HashIntoType>& seen,
                                               bool truncate) const
{
    
    //if (auto graphptr = graph.lock()) {
    if (graph != NULL) {

        // put a 0 on the breadth queue for each element in the starting node queue
        std::queue<unsigned int> breadth_q(std::deque<unsigned int>(node_q.size(), 0));

        unsigned int cur_breadth = 0;
        const unsigned int max_breadth = _tag_density + 1;

        unsigned int total = 0;
        unsigned int nfound = 0;

        KmerFilter filter = [&] (const Kmer& n) -> bool {
            return set_contains(seen, n);
        };
        Traverser traverser(graph, filter);

        while(!node_q.empty()) {

            Kmer node = node_q.front();
            node_q.pop_front();

            unsigned int breadth = breadth_q.front();
            breadth_q.pop();

            // keep track of seen kmers
            seen.insert(node);
            total++;

            // Found a tag!
            if (this->contains(node)) {
                found_tags.push_back(node);
                if (truncate) {
                    return;
                }
                continue;
            }

            if (!(breadth >= cur_breadth)) {
                throw oxli_exception("Desynchonization between traversal "
                                      "and breadth tracking. Did you forget "
                                      "to pop the node or breadth queue?");
            }
            if (breadth > cur_breadth) {
                cur_breadth = breadth;
            }

            if (breadth >= max_breadth) {
                continue;    // truncate search @CTB exit?
            }

            nfound = traverser.traverse(node, node_q);
            for (unsigned int i = 0; i<nfound; ++i) {
                breadth_q.push(breadth + 1);
            }
            total += nfound;
        }
    } else {
        throw oxli_ptr_exception("Hashgraph has been deleted.");
    }
}

