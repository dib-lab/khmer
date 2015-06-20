//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

#include "hashtable.hh"
#include "traversal.hh"

using namespace khmer;
using namespace std;

Traverser::Traverser(const Hashtable * ht) :
    KmerFactory(ht->ksize()), graph(ht)
{
    bitmask = 0;
    for (unsigned int i = 0; i < _ksize; i++) {
        bitmask = (bitmask << 2) | 3;
    }
    rc_left_shift = _ksize * 2 - 2;
}

Kmer Traverser::get_left(Kmer& node, const char ch)
{
    HashIntoType kmer_f, kmer_r;
    kmer_f = ((node.kmer_f) >> 2 | twobit_repr(ch) << rc_left_shift);
    kmer_r = (((node.kmer_r) << 2) & bitmask) | (twobit_comp(ch));
    return build_kmer(kmer_f, kmer_r);
}


Kmer Traverser::get_right(Kmer& node, const char ch)
{
    HashIntoType kmer_f, kmer_r;
    kmer_f = (((node.kmer_f) << 2) & bitmask) | (twobit_repr(ch));
    kmer_r = ((node.kmer_r) >> 2) | (twobit_comp(ch) << rc_left_shift);
    return build_kmer(kmer_f, kmer_r);
}

unsigned int Traverser::traverse_left(Kmer& node,
                                       KmerQueue & target_q,
                                       std::function<bool (Kmer&)> keep_func)
{
    unsigned int found = 0;

    char bases[] = "ACGT";
    char * base = bases;
    while(*base != '\0') {
        Kmer prev_node = get_left(node, *base);
        if (graph->get_count(prev_node) && keep_func(prev_node)) {
            target_q.push_back(prev_node);
            ++found;
        }
        ++base;
    }

    return found;
}

unsigned int Traverser::traverse_right(Kmer& node,
                                       KmerQueue & target_q,
                                       std::function<bool (Kmer&)> keep_func)
{
    unsigned int found = 0;

    char bases[] = "ACGT";
    char * base = bases;
    while(*base != '\0') {
        Kmer next_node = get_right(node, *base);
        if (graph->get_count(next_node) && keep_func(next_node)) {
            target_q.push_back(next_node);
            ++found;
        }
        ++base;
    }

    return found;
}

unsigned int Traverser::degree_left(Kmer& node)
{
    unsigned int degree = 0;

    char bases[] = "ACGT";
    char * base = bases;
    while(*base != '\0') {
        Kmer prev_node = get_left(node, *base);
        if (graph->get_count(prev_node)) {
            ++degree;
        }
        ++base;
    }

    return degree;
}

unsigned int Traverser::degree_right(Kmer& node)
{
    unsigned int degree = 0;

    char bases[] = "ACGT";
    char * base = bases;
    while(*base != '\0') {
        Kmer next_node = get_right(node, *base);
        if (graph->get_count(next_node)) {
            ++degree;
        }
        ++base;
    }

    return degree;
}

unsigned int Traverser::degree(Kmer& node)
{
    return degree_right(node) + degree_left(node);
}

BreadthFirstTraversal::BreadthFirstTraversal(const Hashtable * ht) :
    Traverser(ht)
{
    current_breadth = 0;
    total = 0;
}

unsigned int
BreadthFirstTraversal::search(
       KmerSet& start_nodes,
       KmerSet& start_seen_set,
       std::function<bool ()> continue_func,
       std::function<bool ()> break_func,
       std::function<bool (Kmer& node)> node_keep_func)
{

    auto it = start_nodes.begin();
    Kmer start_node = *it;
    ++it;

    for (; it!=start_nodes.end(); ++it) {
        node_q.push_back(*it);
        breadth_q.push_back(0);
    }

    return search(start_node, start_seen_set, continue_func, break_func,
                  node_keep_func);
}

unsigned int
BreadthFirstTraversal::search(
       Kmer& start_node,
       KmerSet& start_seen_set,
       std::function<bool ()> continue_func,
       std::function<bool ()> break_func,
       std::function<bool (Kmer& node)> node_keep_func)
{
    current_breadth = 0;
    total = 0;
    unsigned int nfound = 0;
    first_node = true;

    seen_set = &start_seen_set;

    node_q.push_front(start_node);
    breadth_q.push_front(0);

    auto filter = [&] (Kmer& n) -> bool {
        return node_keep_func(n);
    };

    int ncontinues = 0;

    while(!node_q.empty()) {

        if (break_func()) {
            break;
        }

        current_node = node_q.front();
        node_q.pop_front();

        current_breadth = breadth_q.front();
        breadth_q.pop_front();

        if (set_contains(*seen_set, current_node)) {
            continue;
        }

        // keep track of seen kmers
        seen_set->insert(current_node);
        total++;

        if (current_breadth > 150
            || total > 10000) {
            std::cout << "Uh oh. Current breadth is " << current_breadth
                      << " and total traversed is " << total
                      << ". Warp core overload imminent!" << std::endl
                      << "Continues: " << ncontinues << std::endl;
            throw khmer_exception("Traversal failed to terminate, "
                                  "warp core has breached. Your ship explodes.");
        }

        if(continue_func()) {
            ++ncontinues;
            continue;
        } else {
            nfound = traverse_right(current_node, node_q, filter);
            for (unsigned int i = 0; i<nfound; ++i) {
                breadth_q.push_back(current_breadth + 1);
            }

            nfound = traverse_left(current_node, node_q, filter);
            for (unsigned int i = 0; i<nfound; ++i) {
                breadth_q.push_back(current_breadth + 1);
            }
        }

        first_node = false;
    }

    /*
    std::cout << std::endl << "Breadth: " << current_breadth << std::endl
              << "Total: " << total << std::endl
              << "Queue size: " << node_q.size() << std::endl
              << "Set size: " << seen_set->size() << std::endl
              << "Breaks: " << nbreaks << std::endl
              << "Continues: " << ncontinues << std::endl;
              */
    return total;
}

bool BreadthFirstTraversal::seen_set_contains(Kmer& node)
{
    return set_contains(*seen_set, node);
}
