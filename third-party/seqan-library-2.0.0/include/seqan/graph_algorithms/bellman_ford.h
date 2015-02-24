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
// Implementation of Bellman-Ford algorithm.
// ==========================================================================

#ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_BELLMAN_FORD_H_
#define INCLUDE_SEQAN_GRAPH_ALGORITHMS_BELLMAN_FORD_H_

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
// Function bellmanFordAlgorithm()
// ----------------------------------------------------------------------------

/*!
 * @fn bellmanFordAlgorithm
 *
 * @headerfile <seqan/graph_algorithms.h>
 *
 * @brief Computes shortest paths from a single source in a directed graph.
 *
 * @signature bool bellmanFordAlgorithm(g, source, weight, predecessor, distance)
 *
 * @param[out] predecessor A property map.  A property map that represents predecessor relationships among vertices.
 *                         It determines a shortest-paths tree.
 * @param[out] distance    A property map.Indicates for each vertex the distance from the source.
 * @param[in]  g           A directed graph. Types: Directed Graph
 * @param[in]  source      A source vertex. Types: VertexDescriptor
 * @param[in]  weight      A weight map.A property map with edge weights.  Edge weights may be negative.
 *
 * @return bool true if the graph has no negative weight cycles, false otherwise.
 *
 * Edge weights may be negative in the Bellman-Ford algorithm.  The out parameters are only valid if the algorithm
 * returns true.
 *
 * @section Example
 *
 * @include demos/graph_algorithms/bellman_ford_algorithm.cpp
 *
 * @include demos/graph_algorithms/bellman_ford_algorithm.cpp.stdout
 *
 * @see dagShortestPath
 * @see dijkstra
 */
template <typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap>
bool bellmanFordAlgorithm(TPredecessorMap & predecessor,
                          TDistanceMap & distance,
                          Graph<TSpec> const & g,
                          TVertexDescriptor const source,
                          TWeightMap const & weight)
{
    typedef typename Size<Graph<TSpec> >::Type TSize;

    // Initialization
    typedef typename Iterator<Graph<TSpec>, VertexIterator>::Type TVertexIterator;
    typedef typename Iterator<Graph<TSpec>, OutEdgeIterator>::Type TOutEdgeIterator;
    resizeVertexMap(predecessor, g);
    resizeVertexMap(distance, g);
    _initializeSingleSource(predecessor, distance, g, source, weight);

    // Run Bellman-Ford
    for(TSize i=0; i<numVertices(g) - 1; ++i) {
        TVertexIterator it(g);
        for(;!atEnd(it);goNext(it)) {
            TVertexDescriptor u = getValue(it);
            TOutEdgeIterator itout(g, u);
            for(;!atEnd(itout);++itout) {
                _relax(g,weight,predecessor, distance, u, getValue(itout));
            }
        }
    }

    TVertexIterator it(g);
    for(;!atEnd(it);goNext(it)) {
        TVertexDescriptor u = getValue(it);
        TOutEdgeIterator itout(g, u);
        for(;!atEnd(itout);++itout) {
            TVertexDescriptor v = targetVertex(g, getValue(itout));
            if (getProperty(distance, v) > getProperty(distance,u) + getProperty(weight,getValue(itout))) {
                return false;
            }
        }
    }
    return true;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_BELLMAN_FORD_H_
