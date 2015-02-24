// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Tobias Raussch <rausch@embl.de>
// ==========================================================================
// Implementation of Depth-First-Search algorithm.
// ==========================================================================

#ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_DEPTH_FIRST_SEARCH_H_
#define INCLUDE_SEQAN_GRAPH_ALGORITHMS_DEPTH_FIRST_SEARCH_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function depthFirstSearch()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TVertexDescriptor, typename TTokenMap, typename TPredecessorMap, typename TDiscoveryTimeMap, typename TFinishingTimeMap, typename TVal>
void
_dfsVisit(Graph<TSpec> const& g,
           TVertexDescriptor const u,
           TTokenMap& tokenMap,
           TPredecessorMap& predecessor,
           TDiscoveryTimeMap& disc,
           TFinishingTimeMap& finish,
           TVal& time)
{
    typedef typename Iterator<Graph<TSpec>, AdjacencyIterator>::Type TAdjacencyIterator;

    assignProperty(tokenMap, u, true);
    ++time;
    assignProperty(disc, u, time);
    TAdjacencyIterator itad(g,u);
    for(;!atEnd(itad);goNext(itad)) {
        TVertexDescriptor v = getValue(itad);
        if (getProperty(tokenMap, v) == false) {
            assignProperty(predecessor, v, u);
            _dfsVisit(g, v, tokenMap, predecessor, disc, finish, time);
        }
    }
    ++time;
    assignProperty(finish, u, time);
}

/*!
 * @fn depthFirstSearch
 * @headerfile <seqan/graph_algorithms.h>
 * @brief Implements a depth-first search on a graph.
 *
 * @signature void depthFirstSearch(predecessor, discovery, finish, g);
 *
 * @param[out] predecessor A property map.Predecessor subgraph produced by the depth-first search.
 * @param[out] discovery   A property map.The discovery time of a vertex v.
 * @param[out] finish      A property map.The time when v's adjacency list has been fully explored.
 * @param[in]  g           A graph. Types: Undirected Graph, Directed Graph
 *
 * In contrast to a breadth-first search the depth-first search is repeated from multiple sources if the graph is not
 * connected.  Hence, depth-first search produces a depth-first forest.  To ensure each vertex ends up in exactly one
 * tree we need not just a distance but a discovery and finishing time.
 *
 * @section Example
 *
 * @include demos/graph_algorithms/depth_first_search.cpp
 *
 * @include demos/graph_algorithms/depth_first_search.cpp.stdout
 *
 * @see breadthFirstSearch
 */
template <typename TSpec, typename TPredecessorMap, typename TDiscoveryTimeMap, typename TFinishingTimeMap>
void depthFirstSearch(TPredecessorMap & predecessor,
                      TDiscoveryTimeMap & disc,
                      TFinishingTimeMap & finish,
                      Graph<TSpec> const & g)
{
    typedef Graph<TSpec> TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Value<TPredecessorMap>::Type TPredVal;

    // Initialization
    resizeVertexMap(predecessor, g);
    resizeVertexMap(disc, g);
    resizeVertexMap(finish, g);
    TPredVal nilPred = getNil<TVertexDescriptor>();

    String<bool> tokenMap;
    resizeVertexMap(tokenMap, g);
    TVertexIterator it(g);
    for(;!atEnd(it);goNext(it)) {
        assignProperty(tokenMap, getValue(it), false);
        assignProperty(predecessor, getValue(it), nilPred);
    }

    TSize time = 0;

    goBegin(it);
    for(;!atEnd(it);goNext(it)) {
        TVertexDescriptor u = getValue(it);
        if (getProperty(tokenMap, u) == false) {
            _dfsVisit(g, u, tokenMap, predecessor, disc, finish, time);
        }
    }
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_DEPTH_FIRST_SEARCH_H_
