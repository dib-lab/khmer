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

Traverser::Traverser(Hashtable * ht) : graph(ht)
{
    bitmask = graph->bitmask;
    K = graph->ksize();
    rc_left_shift = K * 2 - 2;
}

KmerNode Traverser::get_next(KmerNode& node, const char ch)
{
    HashIntoType kmer_f, kmer_r;
    kmer_f = (((node.kmer_f) << 2) & bitmask) | (twobit_repr(ch));
    kmer_r = ((node.kmer_r) >> 2) | (twobit_comp(ch) << rc_left_shift);
    return KmerNode(kmer_f, kmer_r, K);
}

KmerNode Traverser::get_prev(KmerNode& node, const char ch)
{
    HashIntoType kmer_f, kmer_r;
    kmer_f = ((node.kmer_f) >> 2 | twobit_repr(ch) << rc_left_shift);
    kmer_r = (((node.kmer_r) << 2) & bitmask) | (twobit_comp(ch));
    return KmerNode(kmer_f, kmer_r, K);
}

unsigned int Traverser::traverse_right(KmerNode& node,
                                       KmerNodeQueue & node_q)
{
    unsigned int degree = 0;

    char bases[] = "ATCG";
    char * base = bases;
    while(base != NULL) {
        KmerNode next_node = get_next(node, *base);
        if (graph->get_count(next_node.kmer_u)) {
            node_q.push(next_node);
            ++degree;
        }
        ++base;
    }
   
    return degree;
}

unsigned int Traverser::traverse_left(KmerNode& node,
                                       KmerNodeQueue & node_q)
{
    unsigned int degree = 0;

    char bases[] = "ATCG";
    char * base = bases;
    while(base != NULL) {
        KmerNode prev_node = get_prev(node, *base);
        if (graph->get_count(prev_node.kmer_u)) {
            node_q.push(prev_node);
            ++degree;
        }
        ++base;
    }
   
    return degree;
}

unsigned int Traverser::degree_right(KmerNode& node)
{
    unsigned int degree = 0;

    char bases[] = "ATCG";
    char * base = bases;
    while(base != NULL) {
        if (graph->get_count(get_next(node, *base).kmer_u)) {
            ++degree;
        }
        ++base;
    }
   
    return degree;
}

unsigned int Traverser::degree_left(KmerNode& node)
{
    unsigned int degree = 0;

    char bases[] = "ATCG";
    char * base = bases;
    while(base != NULL) {
        if (graph->get_count(get_prev(node, *base).kmer_u)) {
            ++degree;
        }
        ++base;
    }
   
    return degree;
}

unsigned int Traverser::degree(KmerNode& node)
{
    return degree_right(node) + degree_left(node);
}

KmerNode Traverser::build_node(HashIntoType kmer) {
    return KmerNode(kmer, K);
}

KmerNode Traverser::build_node(HashIntoType kmer_f, HashIntoType kmer_r) {
    return KmerNode(kmer_f, kmer_r, K);
}

KmerNode Traverser::build_node(std::string kmer) {
    return KmerNode(kmer, K);
}

KmerNode Traverser::build_node(const char * kmer) {
    return KmerNode(kmer, K);
}
