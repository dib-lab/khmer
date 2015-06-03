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
// Implementation of Kruskal's algorithm.
// ==========================================================================

#ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_KRUSKAL_H_
#define INCLUDE_SEQAN_GRAPH_ALGORITHMS_KRUSKAL_H_

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
// Function kruskalsAlgorithm()
// ----------------------------------------------------------------------------

template <typename TWeight, typename TPair>
struct LessPairI1_ : public std::unary_function<Pair<TWeight, TPair>, bool>
{
    bool operator() (Pair<TWeight, TPair> const & a1,
                     Pair<TWeight, TPair> const & a2) const
    {
        return (a1.i1 < a2.i1);
    }
};

/*!
 * @fn kruskalsAlgorithm
 * @headerfile <seqan/graph_algorithms.h>
 * @brief Computes a minimum spanning tree on a graph.
 *
 * @signature void kruskalsAlgorithm(edges, g, source, weight);
 *
 * @param[out] edges A String of vertex descriptors that represent edges.  Each consecutive pair is an edge with the
 *                   two end points.
 * @param[in] g      An undirected graph. Types: Undirected Graph
 * @param[in] source A source vertex. Types: VertexDescriptor
 * @param[in] weight Edge weights.
 *
 * @section Example
 *
 * @include demos/graph_algorithms/kruskals_algorithm.cpp
 *
 * @include demos/graph_algorithms/kruskals_algorithm.cpp.stdout
 *
 * @see primsAlgorithm
 */
template <typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TEdges>
void kruskalsAlgorithm(TEdges & edges,
                       Graph<TSpec> const & g,
                       TVertexDescriptor const,
                       TWeightMap const & weight)
{
    typedef Graph<TSpec> TGraph;
    typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
    typedef typename Value<TWeightMap>::Type TWeight;

    typedef Pair<TVertexDescriptor, TVertexDescriptor> TVertexPair;
    typedef Pair<TWeight, TVertexPair> TWeightEdgePair;
    typedef String<TWeightEdgePair>  TEdgeList;
    typedef typename Iterator<TEdgeList>::Type TEdgeListIter;
    TEdgeList edgeList;

    // Initialization
    reserve(edges, 2 * (numVertices(g) - 1));
    UnionFind<TVertexDescriptor> unionFind;
    resizeVertexMap(unionFind, g);

    // Sort the edges
    TEdgeIterator itE(g);
    for(;!atEnd(itE);goNext(itE)) appendValue(edgeList, TWeightEdgePair(getProperty(weight, getValue(itE)), TVertexPair(sourceVertex(itE),targetVertex(itE))));
    std::sort(begin(edgeList, Standard() ), end(edgeList, Standard() ), LessPairI1_<TWeight, TVertexPair>() );

    // Process each edge
    TEdgeListIter itEdgeList = begin(edgeList, Standard());
    TEdgeListIter itEdgeListEnd = end(edgeList, Standard());
    for (; itEdgeList != itEdgeListEnd; goNext(itEdgeList)) {
        TVertexDescriptor x = value(itEdgeList).i2.i1;
        TVertexDescriptor y = value(itEdgeList).i2.i2;

        if (findSet(unionFind, x) == findSet(unionFind, y))
            continue;

        appendValue(edges, x);
        appendValue(edges, y);
        joinSets(unionFind, findSet(unionFind, x), findSet(unionFind, y));
    }
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_KRUSKAL_H_
