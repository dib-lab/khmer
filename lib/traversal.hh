/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2015-2016, The Regents of the University of California.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Michigan State University nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
LICENSE (END)

Contact: khmer-project@idyll.org
*/
#ifndef TRAVERSAL_HH
#define TRAVERSAL_HH

#include <queue>
#include <functional>

#include "khmer.hh"
#include "hashtable.hh"
#include "kmer_hash.hh"
#include "kmer_filters.hh"

namespace khmer
{

#ifndef LEFT
#define LEFT 0
#endif
#ifndef RIGHT
#define RIGHT 1
#endif

class Hashtable;
class LabelHash;

class Traverser: public KmerFactory
{
    friend class Hashtable;

protected:

    HashIntoType bitmask;
    unsigned int rc_left_shift;

public:

    const Hashtable * graph;

    explicit Traverser(const Hashtable * ht);

    Kmer get_left(const Kmer& node, const char ch) const;
    Kmer get_right(const Kmer& node, const char ch) const;

    unsigned int traverse_left(Kmer& node,
                               KmerQueue &node_q,
                               KmerFilter filter=0,
                               unsigned short max_neighbors=4) const;
    unsigned int traverse_right(Kmer& node,
                                KmerQueue &node_q,
                                KmerFilter filter=0,
                                unsigned short max_neighbors=4) const;
    unsigned int traverse(Kmer& node,
                          KmerQueue &node_q,
                          KmerFilter filter=0) const {

        unsigned int found;
        found = traverse_left(node, node_q, filter);
        found += traverse_right(node, node_q, filter);
        return found;
    };

    unsigned int degree_left(const Kmer& node) const;
    unsigned int degree_right(const Kmer& node) const;
    unsigned int degree(const Kmer& node) const;
};


template<bool direction>
class AssemblerTraverser: public Traverser
{

protected:

    KmerFilterList filters;

private:

    Kmer get_neighbor(Kmer& node, const char symbol);

public:

    Kmer cursor;

    explicit AssemblerTraverser(const Hashtable * ht,
                             Kmer start_kmer,
                             KmerFilterList filters);

    char next_symbol();
    bool set_cursor(Kmer& node);
    void push_filter(KmerFilter filter);
    KmerFilter pop_filter();
    unsigned int cursor_degree() const;

    std::string join_contigs(std::string& contig_a, std::string& contig_b) const;
};


template<bool direction>
class NonLoopingAT: public AssemblerTraverser<direction>
{
protected:

    SeenSet * visited;

public:

    explicit NonLoopingAT(const Hashtable * ht,
                          Kmer start_kmer,
                          KmerFilterList filters,
                          SeenSet * visited);
    char next_symbol();
};

}
#endif
