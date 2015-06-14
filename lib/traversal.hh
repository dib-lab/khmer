//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

#ifndef TRAVERSAL_HH
#define TRAVERSAL_HH

#include <queue>
#include <functional>

#include "khmer.hh"

#include "khmer_exception.hh"
#include "read_parsers.hh"
#include "kmer_hash.hh"
#include "hashtable.hh"

namespace khmer {

class Kmer;
typedef std::queue<Kmer> KmerQueue;
typedef std::set<Kmer> KmerSet;

class Kmer
{

public:

    HashIntoType kmer_f, kmer_r, kmer_u;

    Kmer(HashIntoType f, HashIntoType r, HashIntoType u)
    {
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

/*
 * At first, I had many contructors defined on Kmer itself. However, it quickly
 * became apparent that this would not work, as the need to pass in a uint8_t
 * for K made the calls ambiguous for. Then I moved it to Traverser, but that
 * really felt like I was shoehorning it in. So here we our, with a factory
 * object. All apologies.
 */
class KmerFactory
{
protected:
    unsigned int K;

public:

    explicit KmerFactory(unsigned int K)
    {
        this->K = K;
    }

    Kmer build_kmer(HashIntoType kmer_u)
    {
        HashIntoType kmer_f, kmer_r;
        _hash(_revhash(kmer_u, K).c_str(), K, kmer_f, kmer_r);
        return Kmer(kmer_f, kmer_r, kmer_u);
    }

    Kmer build_kmer(HashIntoType kmer_f, HashIntoType kmer_r)
    {
        HashIntoType kmer_u = uniqify_rc(kmer_f, kmer_r);
        return Kmer(kmer_f, kmer_r, kmer_u);
    }

    Kmer build_kmer(std::string kmer_s)
    {
        HashIntoType kmer_f, kmer_r, kmer_u;
        kmer_u = _hash(kmer_s.c_str(), K, kmer_f, kmer_r);
        return Kmer(kmer_f, kmer_r, kmer_u);
    }

    Kmer build_kmer(const char * kmer_c)
    {
        HashIntoType kmer_f, kmer_r, kmer_u;
        kmer_u = _hash(kmer_c, K, kmer_f, kmer_r);
        return Kmer(kmer_f, kmer_r, kmer_u);
    }
};

class Traverser: public KmerFactory
{
    friend class Hashtable;

protected:

    HashIntoType bitmask;
    unsigned int rc_left_shift;

public:

    Hashtable * graph;

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
};



};
#endif
