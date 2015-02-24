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
// Implementation of Single-Source-Shortest-Path algorithm.
// ==========================================================================

#ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_SINGLE_SOURCE_SHORTEST_PATH_H_
#define INCLUDE_SEQAN_GRAPH_ALGORITHMS_SINGLE_SOURCE_SHORTEST_PATH_H_

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
// Function singleSourceShortestPath()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TPredecessorMap, typename TVertexDescriptor, typename TNameMap>
inline void
_printPath(Graph<TSpec> const& g,
            TPredecessorMap const& predecessor,
            TVertexDescriptor const source,
            TVertexDescriptor const v,
            TNameMap const& nameMap)
{
    if (source == v) {
        std::cout << getProperty(nameMap, source);
    } else if (getProperty(predecessor, v) == getNil<typename VertexDescriptor<Graph<TSpec> >::Type>()) {
        std::cout << "No path from " << getProperty(nameMap, source) << " to " << getProperty(nameMap, v) << " exists.";
    } else {
        _printPath(g,predecessor, source, getProperty(predecessor, v), nameMap);
        std::cout << "," << getProperty(nameMap, v);
    }
}

template <typename TSpec, typename TPredecessorMap, typename TVertexDescriptor1, typename TVertexDescriptor2>
inline void
_printPath(Graph<TSpec> const& g,
            TPredecessorMap const& predecessor,
            TVertexDescriptor1 const source,
            TVertexDescriptor2 const v)
{
    if (source == v) {
        std::cout << source;
    } else if (getProperty(predecessor, v) == getNil<typename VertexDescriptor<Graph<TSpec> >::Type>()) {
        std::cout << "No path from " << source << " to " << v << " exists.";
    } else {
        _printPath(g,predecessor, source, getProperty(predecessor, v));
        std::cout << "," << v;
    }
}

template <typename TSpec, typename TPredecessorMap, typename TVertexDescriptor1, typename TVertexDescriptor2, typename TEdgeSet>
inline bool
_collectEdges(Graph<TSpec> const& g,
               TPredecessorMap const& predecessor,
               TVertexDescriptor1 const source,
               TVertexDescriptor2 const v,
               TEdgeSet& edgeSet)
{
    if ((TVertexDescriptor1) source == (TVertexDescriptor1) v) {
        return true;
    } else if (getProperty(predecessor, v) == getNil<typename VertexDescriptor<Graph<TSpec> >::Type>()) {
        return false;
    } else {
        edgeSet.insert(findEdge(g, getProperty(predecessor, v), v));
        return _collectEdges(g,predecessor, source, getProperty(predecessor, v), edgeSet);
    }
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec, typename TPredecessorMap, typename TVertexDescriptor, typename TEdgeSet>
inline bool
_collectEdges(Graph<TSpec> const& g,
               TPredecessorMap const& predecessor,
               TVertexDescriptor const source,
               TEdgeSet& edgeSet)
{
    typedef Iterator<Graph<Undirected<> >, VertexIterator>::Type TVertexIterator;
    TVertexIterator it(g);
    for(;!atEnd(it); goNext(it)) {
        if (!_collectEdges(g, predecessor, source, value(it), edgeSet)) {
            edgeSet.clear();
            return false;
        }
    }
    return true;
}

template <typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap>
inline void
_initializeSingleSource(TPredecessorMap & predecessor,
                        TDistanceMap & distance,
                        Graph<TSpec> const & g,
                        TVertexDescriptor const source,
                        TWeightMap const & weight)
{
    typedef Graph<TSpec> TGraph;
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    typedef typename Value<TPredecessorMap>::Type TPredVal;
    typedef typename Value<TWeightMap>::Type TDistVal;
    TPredVal nilPred = getNil<typename VertexDescriptor<TGraph>::Type>();
    TDistVal infDist = _getInfinityDistance(weight);

    TVertexIterator it(g);
    for(;!atEnd(it);goNext(it)) {
        assignProperty(distance, getValue(it), infDist);
        assignProperty(predecessor, getValue(it), nilPred);
    }
    assignProperty(distance, source, 0);
}

template <typename TSpec, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap, typename TVertexDescriptor, typename TEdgeDescriptor>
inline void
_relax(Graph<TSpec> const& g,
        TWeightMap const& weight,
        TPredecessorMap& predecessor,
        TDistanceMap& distance,
        TVertexDescriptor const u,
        TEdgeDescriptor const e)
{
    TVertexDescriptor v = targetVertex(g,e);
    if (getProperty(distance, v) > getProperty(distance,u) + getProperty(weight,e)) {
        assignProperty(distance, v, getProperty(distance,u) + getProperty(weight,e));
        assignProperty(predecessor, v, u);
    }
}

// TODO(holtgrew): This looks broken (distance too high below).

/*!
 * @fn dagShortestPath
 * @headerfile <seqan/graph_algorithms.h>
 * @brief Computes shortest paths from a single source in a directed acyclic graph (DAG).
 *
 * @signature void dagShortestPath(predecessor, distance, g, source, weight);
 *
 * @param[out] predecessor A property map.  A property map that represents predecessor relationships among vertices.
 *                         It determines a shortest-paths tree.
 * @param[out] distance    A property map.  Indicates for each vertex th distance from the source.
 *                         do exist.
 * @param[in]  g           A directed acyclic graph. Types: Directed Graph
 * @param[in]  source      A source vertex. Types: VertexDescriptor
 * @param[in]  weight      A weight map.  In a directed acyclic graph edge weights can be negative because no cycles
 *
 * @section Example
 *
 * @include demos/graph_algorithms/dag_shortest_path.cpp
 *
 * @include demos/graph_algorithms/dag_shortest_path.cpp.stdout
 *
 * @see bellmanFordAlgorithm
 * @see dijkstra
 */
template <typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap>
void dagShortestPath(TPredecessorMap & predecessor,
                     TDistanceMap & distance,
                     Graph<TSpec> const & g,
                     TVertexDescriptor const source,
                     TWeightMap const & weight)
{
    typedef typename Iterator<Graph<TSpec>, OutEdgeIterator>::Type TOutEdgeIterator;
    typedef typename Iterator<String<TVertexDescriptor>, Rooted>::Type TStringIterator;

    // Initialization
    resizeVertexMap(predecessor, g);
    resizeVertexMap(distance, g);

    // Topological sort
    String<TVertexDescriptor> order;
    topologicalSort(order, g);

    _initializeSingleSource(predecessor, distance, g, source, weight);

    //DAG Shortest Paths
    TStringIterator it = begin(order);
    while(!atEnd(it)) {
        TOutEdgeIterator itout(g, getValue(it));
        for(;!atEnd(itout);++itout) {
            _relax(g,weight,predecessor, distance, getValue(it), getValue(itout));
        }
        goNext(it);
    }
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_SINGLE_SOURCE_SHORTEST_PATH_H_
