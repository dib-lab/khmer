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
// Implementation of Prim's algorithm.
// ==========================================================================

#ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_PRIM_H_
#define INCLUDE_SEQAN_GRAPH_ALGORITHMS_PRIM_H_

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
// Function primsAlgorithm()
// ----------------------------------------------------------------------------

/*!
 * @fn primsAlgorithm
 * @headerfile <seqan/graph_algorithms.h>
 * @brief Computes a minimum spanning tree on a graph.
 *
 * @signature void primsAlgorithm(predecessor, g, source, weight);
 *
 * @param[out] predecessor
 *                   A property map.  A property map that represents predecessor relationships among vertices.  It
 *                   determines a minimum spanning tree.
 * @param[in] g      An undirected graph. Types: Undirected Graph
 * @param[in] source A source vertex. Types: VertexDescriptor
 * @param[in] weight Edge weights.
 *
 * @section Example
 *
 * @include demos/graph_algorithms/prims_algorithm.cpp
 *
 * @include demos/graph_algorithms/prims_algorithm.cpp.stdout
 *
 * @see kruskalsAlgorithm
 */
template <typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap>
void primsAlgorithm(TPredecessorMap & predecessor,
                    Graph<TSpec> const & g,
                    TVertexDescriptor const source,
                    TWeightMap const & weight)
{
    typedef Graph<TSpec> TGraph;
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
    typedef typename Value<TPredecessorMap>::Type TPred;
    typedef typename Value<TWeightMap>::Type TWeight;

    typedef std::pair<TWeight, TVertexDescriptor> TWeightVertexPair;
    std::priority_queue<TWeightVertexPair, std::vector<TWeightVertexPair>, std::greater<TWeightVertexPair> > q;

    // Initialization
    String<bool> tokenMap;
    String<TWeight> key;
    TPred nilPred = getNil<typename VertexDescriptor<TGraph>::Type>();
    TWeight infWeight = _getInfinityDistance(weight);
    resizeVertexMap(predecessor, g);
    resizeVertexMap(tokenMap, g);
    resizeVertexMap(key, g);

    TVertexIterator it(g);
    while(!atEnd(it)) {
        TVertexDescriptor u = getValue(it);
        if (u == source) q.push(std::make_pair(0, u));
        assignProperty(predecessor, u, nilPred);
        assignProperty(key, u, infWeight);
        assignProperty(tokenMap, u, false);
        goNext(it);
    }

    assignProperty(key, source, 0);
    while(!q.empty()) {
        TVertexDescriptor u = q.top().second;
        q.pop();
        if (getProperty(tokenMap, u)) continue;
        assignProperty(tokenMap, u, true);
        TOutEdgeIterator itOut(g,u);
        while(!atEnd(itOut)) {
            TVertexDescriptor v = targetVertex(itOut);
            TWeight w = getProperty(weight, getValue(itOut));
            if ((!getProperty(tokenMap, v)) &&
                (w < getProperty(key, v))) {
                    assignProperty(predecessor, v, u);
                    assignProperty(key, v, w);
                    q.push(std::make_pair(w, v));
            }
            goNext(itOut);
        }
    }
}

// TODO(holtgrew): Document.
template <typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap>
void primsAlgorithmSpaceEfficient(TPredecessorMap & predecessor,
                                  Graph<TSpec> const & g,
                                  TVertexDescriptor const source,
                                  TWeightMap const & weight)
{
    typedef Graph<TSpec> TGraph;
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
    typedef typename Value<TPredecessorMap>::Type TPred;
    typedef typename Value<TWeightMap>::Type TWeight;

    // Set-up the priority queue
    typedef Pair<TVertexDescriptor, TWeight> TKeyValue;
    typedef HeapTree<TKeyValue, std::less<TWeight>, KeyedHeap<> > TKeyedHeap;
    TKeyedHeap priorityQueue;

    // Initialization
    String<bool> tokenMap;
    TPred nilPred = getNil<typename VertexDescriptor<TGraph>::Type>();
    TWeight infWeight = _getInfinityDistance(weight);
    resizeVertexMap(predecessor, g);
    resizeVertexMap(tokenMap, g);

    TVertexIterator it(g);
    for(;!atEnd(it);goNext(it)) {
        TVertexDescriptor u = value(it);
        heapInsert(priorityQueue, TKeyValue(u, infWeight));
        assignProperty(predecessor, u, nilPred);
        assignProperty(tokenMap, u, false);
    }
    heapChangeValue(priorityQueue, source, 0);

    // Iterate until queue is empty
    while(!empty(priorityQueue)) {
        TKeyValue kv = heapExtractRoot(priorityQueue);
        TVertexDescriptor u = kv.i1;
        assignProperty(tokenMap, u, true);
        if (kv.i2 == infWeight) continue;
        TOutEdgeIterator itOut(g,u);
        for(;!atEnd(itOut);goNext(itOut)) {
            TVertexDescriptor v = targetVertex(itOut);
            if (getProperty(tokenMap, v)) continue;
            TWeight w = getProperty(weight, getValue(itOut));
            if (w < heapGetValue(priorityQueue, v)) {
                assignProperty(predecessor, v, u);
                heapChangeValue(priorityQueue, v, w);
            }
        }
    }
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_PRIM_H_
