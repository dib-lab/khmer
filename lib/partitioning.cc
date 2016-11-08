#include <queue>
#include <set>
#include <deque>
#include <memory>
#include <functional>

#include "hashtable.hh"
#include "hashbits.hh"
#include "counting.hh"
#include "partitioning.hh"

using namespace khmer;

uint64_t Component::n_created = 0;
bool ComponentPtrCompare::operator() (const ComponentPtr& lhs, 
                                      const ComponentPtr& rhs) const {
    return *lhs < *rhs;
}

inline std::ostream& operator<< (std::ostream& stream, Component& comp) {
    stream << "<Component (id=" << comp.component_id << ", n_tags=" 
           << comp.get_n_tags() << ")>";
    return stream;
}


StreamingPartitioner::StreamingPartitioner(Hashtable * graph)  : 
    graph(graph), _tag_density(DEFAULT_TAG_DENSITY)
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
        tag_component_map = std::unique_ptr<GuardedKmerCompMap>(
                                new GuardedKmerCompMap(graph->ksize(), 
                                                       graph->n_tables(),
                                                       graph_max_table_size / (_tag_density-2)));
        components = std::make_shared<ComponentPtrSet>();
    } else {
        throw khmer_ptr_exception("Hashtable has been deleted.");
    }
}


void StreamingPartitioner::map_tags_to_component(std::set<HashIntoType> tags,
                                                 ComponentPtr& comp)
{
    for (auto tag: tags) {
        tag_component_map->set(tag, comp);
        comp->add_tag(tag);
    }
}


void StreamingPartitioner::consume_sequence(const std::string& seq)
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
    if(graph != NULL) {
        KmerIterator kmers(seq.c_str(), graph->ksize());
        unsigned int since = 1;

        std::set<HashIntoType> tags;
        std::set<HashIntoType> seen;
        std::set<HashIntoType> intersection;
        KmerQueue search_from;

        bool in_known_territory = false;
        bool found_tag_in_territory = false;

        // First check if we overlap any tags
        Kmer kmer = kmers.next();
        tags.insert(kmer); //always tag the first k-mer
        bool is_new_kmer = graph->test_and_set_bits(kmer);
        search_from.push(kmer);

        while(!kmers.done()) {
            bool kmer_tagged = false;

            if (is_new_kmer) {
                // A k-mer from U/G must be searched from for tags,
                // as it could be adjacent to a a k-mer in G/U
                if (in_known_territory && found_tag_in_territory) {
                    // If we had found a tag in the U&G component we just
                    // left, add the component to the seen set.
                    seen.insert(intersection.begin(), intersection.end());
                }
                intersection.clear();

                search_from.push(kmer);
                in_known_territory = false;
                found_tag_in_territory = false;
                ++since;
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

            kmer = kmers.next();
            is_new_kmer = graph->test_and_set_bits(kmer);
        }
        tags.insert(kmer);	// always tag the last k-mer
        search_from.push(kmer);
        intersection.clear();

#if(DEBUG_SP)
        std::cout << "Done iterating k-mers" << std::endl;
        std::cout << tags.size() << " tags in sequence" << std::endl;
#endif
        // Now search for tagged nodes connected to U.
        find_connected_tags(search_from, tags, seen, false);

        // Now resolve components. First, get components from existing tags.
        ComponentPtrSet comps;
        ComponentPtr comp;
#if(DEBUG_SP)
        std::cout << "Get found comps: " << std::endl;
#endif
        for (auto tag: tags) {
#if(DEBUG_SP)
            std::cout << "Tag: " <<  tag;
#endif
            if ((comp = tag_component_map->get(tag)) != NULL) {
#if(DEBUG_SP)
                std::cout << "->" << *comp;
#endif
                comps.insert(comp);
            }
#if(DEBUG_SP)
            std::cout << std::endl;
#endif
        }

#if(DEBUG_SP)
        std::cout << comps.size() << " unique components." << std::endl;
#endif

        if (comps.size() == 0) {
            comp = std::make_shared<Component>();
#if(DEBUG_SP)
            std::cout << "Build new comp: " << *comp << std::endl;
#endif
            components->insert(comp);
        } else {
            // get the first component
            comp = *(comps.begin());
#if(DEBUG_SP)
            std::cout << "Merge into: " << *comp << std::endl;
#endif
                // map the new tags to this component
            comp->add_tag(tags);
            merge_components(comp, comps);
        }
        // (re)map all the tags to the component
        map_tags_to_component(tags, comp);

#if(DEBUG_SP)
        for (auto tag: tags) {
            std::cout << tag << " -> " << *(tag_component_map->get(tag)) << std::endl;
        }
        for (auto c : *components) {
            std::cout << *c << std::endl;
        }
#endif
    } else {
        throw khmer_ptr_exception("Hashtable has been deleted.");
    }
}


void StreamingPartitioner::merge_components(ComponentPtr root, 
                                            ComponentPtrSet comps)
{
#if(DEBUG_SP)
    std::cout << "Merge components." << std::endl;
#endif
    root->merge(comps);
    for (auto other : comps) {
        if (*other == *root) {
            continue;
        }
#if(DEBUG_SP)
        std::cout << "Erase " << *other << std::endl;
#endif
        components->erase(other);
    }
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
        const unsigned int max_breadth = (2 * _tag_density) + 1;

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
                throw khmer_exception("Desynchonization between traversal "
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
        throw khmer_ptr_exception("Hashtable has been deleted.");
    }
}

