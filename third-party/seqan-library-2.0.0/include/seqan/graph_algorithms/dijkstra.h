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
// Implementation of Dijkstra's algorithm.
// ==========================================================================

#ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_DIJKSTRA_H_
#define INCLUDE_SEQAN_GRAPH_ALGORITHMS_DIJKSTRA_H_

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
// Function dijkstra()
// ----------------------------------------------------------------------------

  /*!
 * @fn dijkstra
 * @headerfile <seqan/graph_algorithms.h>
 * @brief Computes shortest paths from a single source in a graph.
 * @signature void dijkstra(predecessor, distance, g, source, weight);
 *
 * @param[out] predecessor A property map.  A property map that represents predecessor relationships among vertices.
 *                         It determines a shortest-paths tree.
 * @param[out] distance    A property map.Indicates for each vertex the distance from the source.
 * @param[in]  g           A graph. Types: Directed Graph
 * @param[in]  source      A source vertex.  Types: VertexDescriptor
 * @param[in]  weight      A weight map.  A property map with edge weights. Edge weights have to be nonnegative.
 *
 * @section Remarks
 *
 * Edge weights have to be nonnegative.
 *
 * @section Example
 *
 * @include demos/graph_algorithms/dijkstra.cpp
 *
 * @include demos/graph_algorithms/dijkstra.cpp.stdout
 *
 * @see dagShortestPath
 * @see bellmanFordAlgorithm
 */
template <typename TSpec, typename TVertexDescriptor, typename TWeightMap,
          typename TPredecessorMap, typename TDistanceMap>
void dijkstra(TPredecessorMap & predecessor,
              TDistanceMap & distance,
              Graph<TSpec> const & g,
              TVertexDescriptor const source,
              TWeightMap const & weight)
{
    typedef Graph<TSpec> TGraph;
    typedef typename Value<TDistanceMap>::Type TDistVal;
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

    // Initialization
    resizeVertexMap(predecessor, g);
    resizeVertexMap(distance, g);

    // S is initially empty
    String<bool> setS;
    resize(setS, getIdUpperBound(_getVertexIdManager(g)), false);

    // Set-up the priority queue
    typedef Pair<TVertexDescriptor, TDistVal> TKeyValue;
    typedef HeapTree<TKeyValue, std::less<TDistVal>, KeyedHeap<> > TKeyedHeap;
    TKeyedHeap priorityQueue;
    TDistVal infDist = _getInfinityDistance(weight);
    TVertexDescriptor nilVertex = getNil<typename VertexDescriptor<TGraph>::Type>();
    TVertexIterator it(g);
    for(;!atEnd(it);goNext(it)) {
        assignProperty(predecessor, value(it), nilVertex);
        assignProperty(distance, value(it), infDist);
        heapInsert(priorityQueue, TKeyValue(value(it), infDist));
    }
    assignProperty(distance, source, 0);
    heapChangeValue(priorityQueue, source, 0);

    // Run Dijkstra
    while (!empty(priorityQueue)) {
        // Extract min
        TVertexDescriptor u = heapExtractRoot(priorityQueue).i1;
        assignProperty(setS, u, true);
        TOutEdgeIterator itout(g, u);
        for(;!atEnd(itout);++itout) {
            TVertexDescriptor v = targetVertex(itout);
            if (property(setS, v) == true) continue;
            if (getProperty(distance, v) > getProperty(distance,u) + getProperty(weight,value(itout))) {
                assignProperty(distance, v, getProperty(distance,u) + getProperty(weight,value(itout)));
                assignProperty(predecessor, v, u);
                heapChangeValue(priorityQueue, v, getProperty(distance,u) + getProperty(weight,value(itout)));
            }
        }
    }
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_DIJKSTRA_H_
