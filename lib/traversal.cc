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
                                       KmerQueue & node_q,
                                       std::function<bool (Kmer&)> filter)
{
    unsigned int found = 0;

    char bases[] = "ACGT";
    char * base = bases;
    while(*base != '\0') {
        Kmer prev_node = get_left(node, *base);
        if (graph->get_count(prev_node) && filter(prev_node)) {
            node_q.push(prev_node);
            ++found;
        }
        ++base;
    }

    return found;
}

unsigned int Traverser::traverse_right(Kmer& node,
                                       KmerQueue & node_q,
                                       std::function<bool (Kmer&)> filter)
{
    unsigned int found = 0;

    char bases[] = "ACGT";
    char * base = bases;
    while(*base != '\0') {
        Kmer next_node = get_right(node, *base);
        if (graph->get_count(next_node) && filter(next_node)) {
            node_q.push(next_node);
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

unsigned int BreadthFirstTraversal::search(Kmer& start_node,
                                           KmerSet& seen_set)
{
    current_breadth = 0;
    total = 0;
    unsigned int nfound = 0;

    this->seen_set = seen_set;

    node_q.push(start_node);
    breadth_q.push(0);

    auto filter = [&] (Kmer& n) -> bool {
        return node_filter_func(n);
    };

    while(!node_q.empty()) {

        if (break_func()) {
            break;
        }

        current_node = node_q.front();
        node_q.pop();

        current_breadth = breadth_q.front();
        breadth_q.pop();

        // keep track of seen kmers
        seen_set.insert(current_node);
        total++;

        if(continue_func()) {
            continue;
        } else {
            nfound = traverse_right(current_node, node_q, filter);
            for (unsigned int i = 0; i<nfound; ++i) {
                breadth_q.push(current_breadth + 1);
            }

            nfound = traverse_left(current_node, node_q, filter);
            for (unsigned int i = 0; i<nfound; ++i) {
                breadth_q.push(current_breadth + 1);
            }
        }

        first_node = false;
    }
}
