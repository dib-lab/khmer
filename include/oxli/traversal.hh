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

#include "oxli.hh"
#include "hashtable.hh"
#include "kmer_hash.hh"
#include "kmer_filters.hh"

namespace oxli
{

#ifndef TRAVERSAL_LEFT
#define TRAVERSAL_LEFT 0
#endif
#ifndef TRAVERSAL_RIGHT
#define TRAVERSAL_RIGHT 1
#endif

class Hashgraph;
class LabelHash;

/**
 * @brief Gather neighbors from a given node.
 *
 * The most basic traversal utility. Stores a list of KmerFilter functions, and given
 * a Kmer, finds all its neighbors that pass the filter function.s
 *
 * @tparam direction The direction in the graph to gather nodes from.
 */
template<bool direction>
class NodeGatherer: public KmerFactory
{
    friend class Hashgraph;

protected:

    KmerFilterList filters;
    HashIntoType bitmask;
    unsigned int rc_left_shift;
    const Hashgraph * graph;

public:

    explicit NodeGatherer(const Hashgraph * ht,
                          KmerFilterList filters);

    explicit NodeGatherer(const Hashgraph * ht);

    explicit NodeGatherer(const Hashgraph * ht, KmerFilter filter);

    WordLength ksize() const;

    /**
     * @brief Push a new filter on to the filter stack.
     */
    void push_filter(KmerFilter filter)
    {
        filters.push_back(filter);
    }

    /**
     * @brief Pop a filter off the stack.
     *
     * @return The filter.
     */
    KmerFilter pop_filter()
    {
        KmerFilter back = this->filters.back();
        this->filters.pop_back();
        return back;
    }

    unsigned int n_filters()
    {
        return filters.size();
    }

    /**
     * @brief Build the Kmer for the potential neighbor node of the given Kmer.
     *
     * When templated for RIGHT, will return the Kmer built from the length K-1 suffix of the
     * input Kmer with the new base appended; when LEFT, the length K-1 prefix of the input Kmer
     * with the new base prepended.
     *
     * @param node The starting node.
     * @param ch The new base to build from.
     *
     * @return The new Kmer.
     */
    Kmer get_neighbor(const Kmer& node, const char ch) const;

    /**
     * @brief Get all neighbors which are present in the graph and pass the filters.
     *
     * @param node The Kmer to start at.
     * @param node_q To collect the results.
     *
     * @return Number of neighbors found.
     */
    unsigned int neighbors(const Kmer& node,
                           KmerQueue &node_q) const;

    /**
     * @brief Get the degree of the given Kmer in the templated direction.
     *
     * @param node The Kmer to check.
     *
     * @return The degree.
     */
    unsigned int degree(const Kmer& node) const;
};


/**
 * @brief A stateful NodeGatherer. Stores its current position.
 *
 * @tparam direction The direction to gather nodes from.
 */
template <bool direction>
class NodeCursor: public NodeGatherer<direction>
{

public:

    // The current position.
    Kmer cursor;
    using NodeGatherer<direction>::push_filter;

    explicit NodeCursor(const Hashgraph * ht,
                        Kmer start_kmer,
                        KmerFilterList filters);

    explicit NodeCursor(const Hashgraph * ht,
                        Kmer start_kmer);

    explicit NodeCursor(const Hashgraph * ht,
                        Kmer start_kmer,
                        KmerFilter filter);

    /**
     * @brief Get the neighbors for the current position.

     *
     * @param node_q To collection the results.
     *
     * @return Number of neighbors found.
     */
    unsigned int neighbors(KmerQueue& node_q) const
    {
        return NodeGatherer<direction>::neighbors(cursor, node_q);
    }

    /**
     * @return Degree of the current cursor position and direction.
     */
    unsigned int cursor_degree() const;

};


/**
 * @brief Wraps a LEFT and RIGHT NodeGatherer.
 */
class Traverser: public KmerFactory
{

protected:

    const Hashgraph * graph;
    NodeGatherer<TRAVERSAL_LEFT> left_gatherer;
    NodeGatherer<TRAVERSAL_RIGHT> right_gatherer;

public:

    explicit Traverser(const Hashgraph * ht,
                       KmerFilterList filters);

    explicit Traverser(const Hashgraph * ht) : Traverser(ht, KmerFilterList()) {}

    explicit Traverser(const Hashgraph * ht,
                       KmerFilter filter);

    void push_filter(KmerFilter filter);
    KmerFilter pop_filter();

    unsigned int traverse(const Kmer& node,
                          KmerQueue& node_q) const;

    unsigned int traverse_left(const Kmer& node,
                               KmerQueue& node_q) const;

    unsigned int traverse_right(const Kmer& node,
                                KmerQueue& node_q) const;

    unsigned int degree(const Kmer& node) const;
    unsigned int degree_left(const Kmer& node) const;
    unsigned int degree_right(const Kmer& node) const;

};


/**
 * @brief A NodeCursor specialized for assembling contigs.
 *
 * @tparam direction The direction to assemble.
 */
template <bool direction>
class AssemblerTraverser: public NodeCursor<direction>
{

protected:
    std::shared_ptr<SeenSet> visited;

public:
    using NodeCursor<direction>::NodeCursor;
    
    explicit AssemblerTraverser(const Hashgraph * ht,
                                Kmer start_kmer,
                                KmerFilterList filters);

    explicit AssemblerTraverser(const Hashgraph * ht,
                                Kmer start_kmer,
                                KmerFilterList filters,
                                std::shared_ptr<SeenSet> visited);

    AssemblerTraverser(const AssemblerTraverser& other);


    /**
     * @brief Get the next symbol.
     *
     * Finds the next symbol which passes the filters, so long as there is only
     * one branch. Does not return a new symbol if there are multiple potential neighbors.
     *
     * @return A member of alphabets::DNA_SIMPLE if a neighbor is found; '\0' otherwise.
     */
    virtual char next_symbol();

    /**
     * @brief Utility function to join two overlapping contigs with proper directionality.
     *
     *  By default, assumes the two contigs overlap by length K. This can be reduced via the
     *  offset parameter.
     *
     * @param contig_a
     * @param contig_b
     * @param offset Number of bases to subtract from K when joining.
     *
     * @return The joined contig.
     */
    std::string join_contigs(std::string& contig_a,
                             std::string& contig_b,
                             WordLength offset = 0) const;
};


}
#endif
