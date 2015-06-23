//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

#ifndef TRAVERSAL_HH
#define TRAVERSAL_HH

#include <deque>
#include <queue>
#include <functional>

#include "khmer.hh"

#include "khmer_exception.hh"
#include "read_parsers.hh"
#include "kmer_hash.hh"
#include "hashtable.hh"

namespace khmer {

class Hashtable;

class Traverser: public KmerFactory
{
    friend class Hashtable;

protected:

    HashIntoType bitmask;
    unsigned int rc_left_shift;

public:

    const Hashtable * graph;

    explicit Traverser(const Hashtable * ht);

    Kmer get_left(Kmer& node, const char ch);
    Kmer get_right(Kmer& node, const char ch);

    unsigned int traverse_left(Kmer& node,
                               KmerQueue &node_q,
                               std::function<bool (Kmer&)> keep_func);
    unsigned int traverse_right(Kmer& node,
                                KmerQueue &node_q,
                                std::function<bool (Kmer&)> keep_func);

    unsigned int degree_left(Kmer& node);
    unsigned int degree_right(Kmer& node);
    unsigned int degree(Kmer& node);
};

class BreadthFirstTraversal: public Traverser
{
protected:

    unsigned int current_breadth, total;
    bool first_node;
    Kmer current_node;

public:

    KmerSet * seen_set;
    KmerQueue node_q;
    std::deque<unsigned int> breadth_q;

    explicit BreadthFirstTraversal(const Hashtable * ht);

    unsigned int
    search(Kmer& start_node,
           KmerSet& start_seen_set,
           std::function<bool ()> continue_func,
           std::function<bool ()> break_func,
           std::function<bool (Kmer& node)> node_keep_func);

   unsigned int
   search(KmerQueue& start_nodes,
          KmerSet& start_seen_set,
          std::function<bool ()> continue_func,
          std::function<bool ()> break_func,
          std::function<bool (Kmer& node)> node_keep_func);

    Kmer cursor() {
        return current_node;
    }

    bool on_first_node() {
        return first_node;
    }

    unsigned int get_cursor_breadth() {
        return current_breadth;
    }

    unsigned int get_num_traversed() {
        return total;
    }

    bool seen_set_contains(Kmer& node);
    unsigned int seen_set_size() {
        return seen_set->size();
    }

};

};
#endif
