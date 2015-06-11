//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

#ifndef TRAVERSAL_HH
#define TRAVERSAL_HH

#include <queue>

#include "khmer.hh"

#include "khmer_exception.hh"
#include "read_parsers.hh"
#include "kmer_hash.hh"
#include "hashtable.hh"

namespace khmer {

class KmerNode;
typedef std::queue<KmerNode> KmerNodeQueue;

class KmerNode {

public:

    HashIntoType kmer_f, kmer_r, kmer_u;

    KmerNode(HashIntoType f, HashIntoType r, HashIntoType u) {
        kmer_f = f;
        kmer_r = r;
        kmer_u = u;
    }

    std::string get_string_rep(unsigned int K) {
        return _revhash(kmer_u, K);
    }

    const char * get_char_rep(unsigned int K) {
        return _revhash(kmer_u, K).c_str();
    }
};

class Traverser {

    friend class Hashtable;
    
protected:

    HashIntoType bitmask;
    unsigned int rc_left_shift;

public:

    Hashtable * graph;
    unsigned int K;

    explicit Traverser(Hashtable * ht);

    KmerNode get_next(KmerNode& node, const char ch);
    KmerNode get_prev(KmerNode& node, const char ch);

    unsigned int traverse_right(KmerNode& node,
                                KmerNodeQueue &node_q);
    unsigned int traverse_left(KmerNode& node,
                               KmerNodeQueue &node_q);
    
    unsigned int degree_right(KmerNode& node);
    unsigned int degree_left(KmerNode& node);
    unsigned int degree(KmerNode& node);

    KmerNode build_node(HashIntoType kmer);
    KmerNode build_node(HashIntoType kmer_f, HashIntoType kmer_r);
    KmerNode build_node(std::string kmer);
    KmerNode build_node(const char * kmer);

};
};
#endif
