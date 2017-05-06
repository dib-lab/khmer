#include <queue>
#include <set>
#include <deque>
#include <memory>
#include <functional>
#include <algorithm>

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


StreamingPartitioner::StreamingPartitioner(Hashgraph * graph, uint32_t tag_density)  : 
    graph(graph), _tag_density(tag_density), components_lock(0), n_consumed(0)
{

    //if (auto graphptr = graph.lock()) {
    if (graph != NULL) {
        std::vector<uint64_t> graph_table_sizes = graph->get_tablesizes(); 
        uint64_t graph_max_table_size = *std::max_element(graph_table_sizes.begin(),
                                                          graph_table_sizes.end());

        // We can guess that, given N k-mers in the graph, there will be
        // approximately N / _tag_density tags. If we want to the filter false
        // positive rate to be about the same as the graph, we should make its table
        // sizes proportional by the number of tags. Here, we use _tag_density-2
        // because we always tag the first and last k-mers in a read.
        tag_component_map = std::unique_ptr<GuardedHashCompMap>(
                                new GuardedHashCompMap(graph->ksize(), 
                                                       graph->n_tables(),
                                                       graph_max_table_size / (_tag_density-2)));
        components = std::make_shared<ComponentPtrSet>();
    } else {
        throw oxli_ptr_exception("Hashgraph has been deleted.");
    }
}


void StreamingPartitioner::map_tags_to_component(std::set<HashIntoType>& tags,
                                                 ComponentPtr& comp)
{
    for (auto tag: tags) {
        tag_component_map->set(tag, comp);
        comp->add_tag(tag);
    }
}


uint64_t StreamingPartitioner::consume_fasta(const std::string& filename)
{
    ReadParserPtr<FastxReader> parser = get_parser<FastxReader>(filename);
    Read read;
    uint64_t n_invalid = 0;
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

        bool is_valid = graph->check_and_normalize_read(read.sequence);
        if (is_valid) {
            consume(read.sequence);
        } else {
            n_invalid++;
        }
        n_consumed++;
    }

    return n_invalid;
}



uint64_t StreamingPartitioner::consume(const std::string& seq)
{
    std::set<HashIntoType> tags;
    KmerQueue seeds;
    std::set<HashIntoType> seen;

    uint64_t n_new = seed_sequence(seq, tags, seeds, seen);
    find_connected_tags(seeds, tags, seen, false);
    //acquire_components();
    create_and_connect_components(tags);
    //release_components();
    return n_new;
}


uint64_t StreamingPartitioner::consume_pair(const std::string& first,
                                        const std::string& second)
{
    std::set<HashIntoType> tags;
    KmerQueue seeds;
    std::set<HashIntoType> seen;

    uint64_t n_new = seed_sequence(first, tags, seeds, seen);
    n_new += seed_sequence(second, tags, seeds, seen);
    find_connected_tags(seeds, tags, seen, false);
    //acquire_components();
    create_and_connect_components(tags);
    //release_components();
    return n_new;
}

void StreamingPartitioner::add_component(ComponentPtr comp)
{
    components->insert(comp);
    map_tags_to_component(comp->tags, comp);
}



uint64_t StreamingPartitioner::seed_sequence(const std::string& seq,
                                          std::set<HashIntoType>& tags,
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

                seeds.push(kmer);
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
                kmer_tagged = tag_component_map->contains(kmer);
                if (kmer_tagged) {
                    since = 1;
                    tags.insert(kmer);
                    found_tag_in_territory = true;
                } else {
                    ++since;
                }
            }

            if (since >= _tag_density) {
                tags.insert(kmer);
                since = 1;
            }
        } while (!kmers.done());

        // always tag the last k-mer
        if (since >= _tag_density / 2) {
            tags.insert(kmer);
        }
        seeds.push(kmer);

        // now go back and make sure to search from the first k-mer
        kmer = kmers.first();
        seeds.push(kmer);

#if(DEBUG_SP)
        std::cout << "Done iterating k-mers" << std::endl;
        std::cout << tags.size() << " tags in sequence" << std::endl;
#endif
    } else {
        throw oxli_ptr_exception("Hashgraph has been deleted.");
    }

    return n_new;
}


uint32_t StreamingPartitioner::create_and_connect_components(std::set<HashIntoType> &tags)
{

        // Now resolve components. First, get components from existing tags.
        ComponentPtrSet found_comps;
#if(DEBUG_SP)
        std::cout << "Get found comps: " << std::endl;
#endif
        for (auto tag: tags) {
#if(DEBUG_SP)
            std::cout << "Tag: " <<  tag;
#endif
            ComponentPtr comp;
            if ((comp = tag_component_map->get(tag)) != NULL) {
#if(DEBUG_SP)
                std::cout << "->" << *comp;
#endif
                found_comps.insert(comp);
            }
#if(DEBUG_SP)
            std::cout << std::endl;
#endif
        }

#if(DEBUG_SP)
        std::cout << found_comps.size() << " unique components." << std::endl;
#endif
        
        uint32_t n_merged = 1;
        if (found_comps.size() == 0) {
            ComponentPtr new_comp = std::make_shared<Component>();
#if(DEBUG_SP)
            std::cout << "Build new comp: " << *new_comp << std::endl;
#endif
            components->insert(new_comp);
            map_tags_to_component(tags, new_comp);
        } else {
            // Choose the largest component as the root
            // We want to minimize tag copying
            ComponentPtr root_comp = *(found_comps.begin());
            for (auto other : found_comps) {
                if (other->get_n_tags() > root_comp->get_n_tags()) {
                    root_comp = other;
                }
            }
#if(DEBUG_SP)
            std::cout << "Merge into: " << *root_comp << std::endl;
#endif
            // map the new tags to this component
            root_comp->add_tag(tags);
            map_tags_to_component(tags, root_comp);
            if (found_comps.size() > 1) {
                n_merged = merge_components(root_comp, found_comps);
            }
        }
        return n_merged;
}


uint32_t StreamingPartitioner::merge_components(ComponentPtr& root, 
                                            ComponentPtrSet& comps)
{
    uint32_t n_merged = 1;
    for (auto other : comps) {
        if (*other == *root) {
            continue;
        }
        root->add_tag(other->tags); // transfer the tags from the other comp
        map_tags_to_component(other->tags, root); // set the other's tags to point to root
        components->erase(other); // remove other component entirely
        n_merged++;

    }
    comps.clear(); // should call destructor on all the merged comps, unless they have
                   // and active Python wrapper; this leaves them as sole owners
    return n_merged;
}


ComponentPtr StreamingPartitioner::get_tag_component(HashIntoType tag) const
{
    return tag_component_map->get(tag);
}


ComponentPtr StreamingPartitioner::get_tag_component(std::string& kmer) const
{
    HashIntoType h = graph->hash_dna(kmer.c_str());
    return get_tag_component(h);
}


ComponentPtr StreamingPartitioner::get_nearest_component(std::string& kmer) const
{
    Kmer hashed = graph->build_kmer(kmer);
    return get_nearest_component(hashed);
}


ComponentPtr StreamingPartitioner::get_nearest_component(Kmer kmer) const
{
    std::set<HashIntoType> tags;
    std::set<HashIntoType> seen;
    KmerQueue node_q;
    node_q.push(kmer);

    find_connected_tags(node_q, tags, seen, true);
    if (tags.size() > 0) {
        HashIntoType tag = *(tags.begin());
        return tag_component_map->get(tag);
    } else {
        return NULL;
    }
}


void StreamingPartitioner::find_connected_tags(KmerQueue& node_q,
                                               std::set<HashIntoType>& found_tags,
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
            node_q.pop();

            unsigned int breadth = breadth_q.front();
            breadth_q.pop();

            // keep track of seen kmers
            seen.insert(node);
            total++;

            // Found a tag!
            if (tag_component_map->contains(node)) {
                found_tags.insert(node);
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

