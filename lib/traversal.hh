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
