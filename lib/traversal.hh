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

class Kmer;
typedef std::queue<Kmer> KmerQueue;
typedef std::set<Kmer> KmerSet;

class Kmer {

public:

    HashIntoType kmer_f, kmer_r, kmer_u;

    Kmer(HashIntoType f, HashIntoType r, HashIntoType u) {
        kmer_f = f;
        kmer_r = r;
        kmer_u = u;
    }

    bool operator< (const Kmer &other) const {
        return kmer_u < other.kmer_u;
    }

    std::string get_string_rep(unsigned int K) {
        return _revhash(kmer_u, K);
    }

    const char * get_char_rep(unsigned int K) {
        return _revhash(kmer_u, K).c_str();
    }
};


class Traverser
{
    friend class Hashtable;

protected:

    HashIntoType bitmask;
    unsigned int rc_left_shift;

public:

    Hashtable * graph;
    unsigned int K;

    explicit Traverser(Hashtable * ht);

    Kmer get_left(Kmer& node, const char ch);
    Kmer get_right(Kmer& node, const char ch);

    unsigned int traverse_left(Kmer& node,
                               KmerQueue &node_q,
                               std::function<bool (Kmer&)> filter);
    unsigned int traverse_right(Kmer& node,
                                KmerQueue &node_q,
                                std::function<bool (Kmer&)> filter);

    unsigned int degree_left(Kmer& node);
    unsigned int degree_right(Kmer& node);
    unsigned int degree(Kmer& node);

    Kmer build_node(HashIntoType kmer);
    Kmer build_node(HashIntoType kmer_f, HashIntoType kmer_r);
    Kmer build_node(std::string kmer);
    Kmer build_node(const char * kmer);

};



};
#endif
