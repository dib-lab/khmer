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
// Implementation of Breadth-First-Search.
// ==========================================================================

#ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_BREADTH_FIRST_SEARCH_H_
#define INCLUDE_SEQAN_GRAPH_ALGORITHMS_BREADTH_FIRST_SEARCH_H_

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
// Function breadthFirstSearch()
// ----------------------------------------------------------------------------

/*!
 * @fn breadthFirstSearch
 * @headerfile <seqan/graph_algorithms.h>
 * @brief Implements a breadth-first search on a graph.
 *
 * @signature void breadthFirstSearch(predecessor, distance, g, source);
 *
 * @param[out] predecessor A property map.  The predecessor map stores implicitly the breadth-first tree.
 * @param[out] distance    A property map.  The distance map indicates at what depth a vertex was discovered.
 * @param[in]  g           Undirected Graph, Directed Graph
 * @param[in]  source      A vertex descriptor.  The breadth-first search is started from this vertex.
 *                         Types: VertexDescriptor
 *
 * Breadth-first search computes the distance from source to all reachable vertices.  It also produces a breath-first
 * tree where each node has a predecessor/parent.
 *
 * @section Example
 *
 * @include demos/graph_algorithms/breadth_first_search.cpp
 *
 * @include demos/graph_algorithms/breadth_first_search.cpp.stdout
 *
 * @see depthFirstSearch
 */

template <typename TSpec, typename TVertexDescriptor, typename TPredecessorMap, typename TDistanceMap>
void breadthFirstSearch(TPredecessorMap & predecessor,
                        TDistanceMap & distance,
                        Graph<TSpec> const & g,
                        TVertexDescriptor const source)
{
    typedef Graph<TSpec> TGraph;
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    typedef typename Value<TPredecessorMap>::Type TPredVal;
    typedef typename Value<TDistanceMap>::Type TDistVal;

    // Initialization
    resizeVertexMap(predecessor, g);
    resizeVertexMap(distance, g);
    TPredVal nilPred = getNil<typename VertexDescriptor<TGraph>::Type>();
    TDistVal infDist = _getInfinityDistance(distance);

    String<bool> tokenMap;
    resizeVertexMap(tokenMap, g);
    TVertexIterator it(g);
    for(;!atEnd(it);goNext(it)) {
        assignProperty(tokenMap, getValue(it), false);
        assignProperty(distance, getValue(it), infDist);
        assignProperty(predecessor, getValue(it), nilPred);
    }
    assignProperty(tokenMap, source, true);
    assignProperty(distance, source, 0);
    assignProperty(predecessor, source, nilPred);
    std::deque<TVertexDescriptor> queue;
    queue.push_back(source);

    // Bfs
    while (!queue.empty()) {
        TVertexDescriptor u = queue.front();
        queue.pop_front();
        typedef typename Iterator<Graph<TSpec>, OutEdgeIterator>::Type TOutEdgeIterator;
        TOutEdgeIterator itout(g,u);
        for(;!atEnd(itout);goNext(itout)) {
            TVertexDescriptor v = targetVertex(itout);
            if (getProperty(tokenMap, v) == false) {
                assignProperty(tokenMap, v, true);
                assignProperty(distance, v, getProperty(distance,u) + 1);
                assignProperty(predecessor, v, u);
                queue.push_back(v);
            }
        }
    }
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_BREADTH_FIRST_SEARCH_H_
