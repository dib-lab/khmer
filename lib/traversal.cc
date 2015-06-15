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

Traverser::Traverser(Hashtable * ht) : KmerFactory(ht->ksize()), graph(ht)
{
    bitmask = graph->bitmask;
    rc_left_shift = ksize * 2 - 2;
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

    char bases[] = "ATCG";
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

    char bases[] = "ATCG";
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

    char bases[] = "ATCG";
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

    char bases[] = "ATCG";
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
