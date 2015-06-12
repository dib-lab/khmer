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
    return build_node(kmer_f, kmer_r);
}

KmerNode Traverser::get_prev(KmerNode& node, const char ch)
{
    HashIntoType kmer_f, kmer_r;
    kmer_f = ((node.kmer_f) >> 2 | twobit_repr(ch) << rc_left_shift);
    kmer_r = (((node.kmer_r) << 2) & bitmask) | (twobit_comp(ch));
    return build_node(kmer_f, kmer_r);
}

unsigned int Traverser::traverse_right(KmerNode& node,
                                       KmerNodeQueue & node_q)
{
    unsigned int degree = 0;

    char bases[] = "ATCG";
    char * base = bases;
    while(*base != '\0') {
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
    while(*base != '\0') {
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
    while(*base != '\0') {
        KmerNode next = get_next(node, *base);
        if (graph->get_count(next.kmer_u)) {
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
    while(*base != '\0') {
        KmerNode prev = get_prev(node, *base);
        if (graph->get_count(prev.kmer_u)) {
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

KmerNode Traverser::build_node(HashIntoType kmer_u) {
    HashIntoType kmer_f, kmer_r;
    _hash(_revhash(kmer_u, K).c_str(), K, kmer_f, kmer_r);
    return KmerNode(kmer_f, kmer_r, kmer_u);
}

KmerNode Traverser::build_node(HashIntoType kmer_f, HashIntoType kmer_r) {
    HashIntoType kmer_u = uniqify_rc(kmer_f, kmer_r);
    return KmerNode(kmer_f, kmer_r, kmer_u);
}

KmerNode Traverser::build_node(std::string kmer) {
    HashIntoType kmer_f, kmer_r, kmer_u;
    kmer_u = _hash(kmer.c_str(), K, kmer_f, kmer_r);
    return KmerNode(kmer_f, kmer_r, kmer_u);
}

KmerNode Traverser::build_node(const char * kmer) {
    HashIntoType kmer_f, kmer_r, kmer_u;
    kmer_u = _hash(kmer, K, kmer_f, kmer_r);
    return KmerNode(kmer_f, kmer_r, kmer_u);
}
