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
// Implementation of the Ford-Fulkerson algorithm.
// ==========================================================================

#ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_FORD_FULKERSON_H_
#define INCLUDE_SEQAN_GRAPH_ALGORITHMS_FORD_FULKERSON_H_

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
// Function fordFulkersonAlgorithm()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TCapMap, typename TFlowMap, typename TResidualGraph>
void
_buildResidualGraph(Graph<TSpec> const& g,
                      TCapMap const& capacity,
                      TFlowMap const& flow,
                      TResidualGraph& rG)
{
    typedef Graph<TSpec> TGraph;
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
    typedef typename Value<TFlowMap>::Type TFlow;
    typedef typename Value<TCapMap>::Type TCap;

    clear(rG);
    TVertexIterator itV(g);
    for(;!atEnd(itV);goNext(itV)) {
        _createVertices(rG, getValue(itV));
    }

    TEdgeIterator itE(g);
    for(;!atEnd(itE);goNext(itE)) {
        typedef typename EdgeDescriptor<TResidualGraph>::Type TEdgeDescriptor;
        TFlow f = getProperty(flow, getValue(itE));
        TCap cap = getProperty(capacity, getValue(itE));
        if (f > 0) {
            TEdgeDescriptor e_rG = findEdge(rG, targetVertex(itE), sourceVertex(itE));
            if (e_rG == 0) addEdge(rG, targetVertex(itE), sourceVertex(itE), f);
            else cargo(e_rG) += f;
        }
        if (f < cap) {
            TEdgeDescriptor e_rG = findEdge(rG, sourceVertex(itE), targetVertex(itE));
            if (e_rG == 0) addEdge(rG, sourceVertex(itE), targetVertex(itE), cap - f);
            else cargo(e_rG) += cap - f;
        }
    }
}

template <typename TSpec, typename TPredecessorMap, typename TVertexDescriptor>
inline typename Size<Graph<TSpec> >::Type
_getMinimumAug(Graph<TSpec> const & rG,
                 TPredecessorMap & predecessor,
                 TVertexDescriptor const source,
                 TVertexDescriptor sink)
{
    typedef Graph<TSpec> TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef TSize TFlow;
    typedef typename Iterator<String<TVertexDescriptor>, Rooted>::Type TIterator;

    // Build secondary predecessor map just containing the path
    TVertexDescriptor nilPred = getNil<typename VertexDescriptor<TGraph>::Type>();
    String<TVertexDescriptor> predMap;
    resizeVertexMap(predMap, rG);
    TIterator it = begin(predMap);
    for(;!atEnd(it);goNext(it)) {
        *it = nilPred;
    }

    // Find minimum flow
    TVertexDescriptor pred = getProperty(predecessor, sink);
    TFlow f = getCargo(findEdge(rG, pred,sink));
    assignProperty(predMap, sink, pred);
    while(pred != source) {
        sink = pred;
        pred = getProperty(predecessor, sink);
        TFlow f2 = getCargo(findEdge(rG, pred,sink));
        assignProperty(predMap, sink, pred);
        if (f2 < f) f = f2;
    }

    // Just return the augmenting path
    predecessor = predMap;
    return f;
}

/*!
 * @fn fordFulkersonAlgorithm
 *
 * @headerfile <seqan/graph_algorithms.h>
 *
 * @brief Computes a maximum flow in a directed graph.
 *
 * @signature TValue fordFulkeronAlgorithm(flow, g, source, sink, capacity);
 *
 * @param[out] flow     A property map with the flow of each edge.
 * @param[in]  g        A directed graph.  Types: Directed Graph
 * @param[in]  source   A source vertex.  Types: VertexDescriptor
 * @param[in]  sink     A sink vertex.  Types: VertexDescriptor
 * @param[in]  capacity A property map of edge capacities.
 *
 * @return TValue The value of the flow.  TValue is the Value tpye of the type of flow.
 *
 * @section Example
 *
 * @include demos/graph_algorithms/ford_fulkerson_algorithm.cpp
 *
 * @include demos/graph_algorithms/ford_fulkerson_algorithm.cpp.stdout
 */
template <typename TSpec, typename TVertexDescriptor, typename TCapMap, typename TFlowMap>
typename Value<TFlowMap>::Type
fordFulkersonAlgorithm(TFlowMap & flow,
                       Graph<TSpec> const & g,
                       TVertexDescriptor const source,
                       TVertexDescriptor const sink,
                       TCapMap const & capacity)
{
    typedef Graph<TSpec> TGraph;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
    typedef typename Value<TFlowMap>::Type TFlow;

    // Initialization
    TVertexDescriptor nilPred = getNil<typename VertexDescriptor<TGraph>::Type>();
    resizeEdgeMap(flow, g);
    TEdgeIterator itE(g);
    for(;!atEnd(itE);goNext(itE)) {
        assignProperty(flow, getValue(itE), 0);
    }

    // Build the residual graph
    Graph<Directed<TFlow> > rG;
    _buildResidualGraph(g,capacity, flow, rG);


    // Determine whether the sink is reachable
    String<TVertexDescriptor> predMap;
    String<TVertexDescriptor> distMap;
    breadthFirstSearch(predMap, distMap, rG, source);

    while (getProperty(predMap, sink) != nilPred) {
        TFlow inc = _getMinimumAug(rG, predMap, source, sink);
        TEdgeIterator itEdge(g);
        for(;!atEnd(itEdge);goNext(itEdge)) {
            TVertexDescriptor u = sourceVertex(itEdge);
            TVertexDescriptor v = targetVertex(itEdge);
            TEdgeDescriptor e = getValue(itEdge);
            if (getProperty(predMap, v) == u) assignProperty(flow, e, getProperty(flow, e) + inc);
            if (getProperty(predMap, u) == v) assignProperty(flow, e, getProperty(flow, e) - inc);
        }
        // Build the residual graph
        _buildResidualGraph(g, capacity, flow, rG);
        // Determine whether the sink is reachable
        clear(predMap);
        clear(distMap);
        breadthFirstSearch(predMap, distMap, rG, source);
    }

    TFlow valF = 0;
    TOutEdgeIterator itOutEdge(g, source);
    for(;!atEnd(itOutEdge);goNext(itOutEdge)) {
        valF += getProperty(flow, getValue(itOutEdge));
    }
    return valF;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_FORD_FULKERSON_H_
