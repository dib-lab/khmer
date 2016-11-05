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

StreamingPartitioner::StreamingPartitioner(Hashtable * graph)  : 
    graph(graph), _tag_density(DEFAULT_TAG_DENSITY), n_components(0)
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
    } else {
        throw khmer_ptr_exception("Hashtable has been deleted.");
    }
}


void StreamingPartitioner::map_tags_to_component(std::set<HashIntoType> tags,
                                                 ComponentPtr& comp)
{
    for (auto tag: tags) {
        tag_component_map->set(tag, comp);
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

            is_new_kmer = graph->test_and_set_bits(kmer);
            kmer = kmers.next();
        }
        tags.insert(kmer);	// always tag the last k-mer
        if (is_new_kmer || !found_tag_in_territory) {
            search_from.push(kmer);
        }
        intersection.clear();

        // Now search for tagged nodes connected to U.
        find_connected_tags(search_from, tags, seen);

        // Now resolve components. First, get components from existing tags.
        std::set<ComponentPtr> comps;
        ComponentPtr comp;
        for (auto tag: tags) {
            if ((comp = tag_component_map->get(tag)) != NULL) {
                comps.insert(comp);
            }
        }

        if (comps.size() == 0) {
            comp = std::make_shared<Component>();
            n_components++;
        } else {
            // get the first component
            comp = *(comps.begin());
            if (comps.size() == 1) {
                // map the new tags to this component
                comp->add_tag(tags);
                map_tags_to_component(tags, comp);
            } else {
                // merge the components
                comp->merge(comps);
                n_components -= comps.size() - 1;
            }
        }
        // (re)map all the tags to the component
        map_tags_to_component(tags, comp);
    } else {
        throw khmer_ptr_exception("Hashtable has been deleted.");
    }
}


void StreamingPartitioner::find_connected_tags(KmerQueue& node_q,
                                               std::set<HashIntoType>& found_tags,
                                               std::set<HashIntoType>& seen)
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
        }
    } else {
        throw khmer_ptr_exception("Hashtable has been deleted.");
    }
}

