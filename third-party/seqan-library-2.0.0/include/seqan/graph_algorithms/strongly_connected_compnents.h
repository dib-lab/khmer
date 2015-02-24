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
// Implementation of Strongly-Connected-Compnonents algorithm.
// ==========================================================================

#ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_STRONGLY_CONNECTED_COMPNENTS_H_
#define INCLUDE_SEQAN_GRAPH_ALGORITHMS_STRONGLY_CONNECTED_COMPNENTS_H_

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
// Function stronglyConnectedComponents()
// ----------------------------------------------------------------------------

/*!
 * @fn stronglyConnectedComponents
 * @headerfile <seqan/graph_algorithms.h>
 * @brief Decomposes a directed graph into its strongly connected components.
 * @signature TSize stronglyConnectedComponents(components, g);
 *
 * @param[out] components
 *               A property map.Each vertex is mapped to a component id.  If two vertices share the same id they are
 *               in the same component.
 * @param[in]  g A directed graph. Types: Directed Graph
 *
 * @return TSize Number of strongly connected components, Size type of g.
 *
 * @section Example
 *
 * @include demos/graph_algorithms/strongly_connected_components.cpp
 *
 * @include demos/graph_algorithms/strongly_connected_components.cpp.stdout
 */
template <typename TSpec, typename TComponents>
typename Size<Graph<TSpec> >::Type
stronglyConnectedComponents(TComponents & components,
                            Graph<TSpec> const & g_source)
{
    // Initialization
    typedef Graph<TSpec> TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    typedef typename Value<TComponents>::Type TCompVal;
    resizeVertexMap(components, g_source);
    String<TSize> predMap;
    String<TSize> discoveryTimeMap;
    String<TSize> finishingTimeMap;

    // Dfs
    depthFirstSearch(predMap, discoveryTimeMap, finishingTimeMap, g_source);

    Graph<TSpec> g;
    transpose(g_source, g);

    // Second Dfs
    String<TSize> predecessor;
    String<TSize> disc;
    String<TSize> finish;
    resizeVertexMap(predecessor, g);
    resizeVertexMap(disc, g);
    resizeVertexMap(finish, g);
    TCompVal nilPred = getNil<TVertexDescriptor>();
    String<bool> tokenMap;
    resizeVertexMap(tokenMap, g);
    TVertexIterator it(g);
    for(;!atEnd(it);goNext(it)) {
        assignProperty(components, getValue(it), nilPred);
        assignProperty(tokenMap, getValue(it), false);
        assignProperty(predecessor, getValue(it), nilPred);
    }

    // Order vertices
    typedef std::pair<TSize, TVertexDescriptor> TTimeVertexPair;
    std::priority_queue<TTimeVertexPair> q;
    goBegin(it);
    for(;!atEnd(it);++it) {
        q.push(std::make_pair(getProperty(finishingTimeMap, getValue(it)), getValue(it)));
    }

    TSize time = 0;
    TSize label = 0;
    while(!q.empty()) {
        TVertexDescriptor u = q.top().second;
        q.pop();
        if (getProperty(tokenMap, u) == false) {
            _dfsVisit(g, u, tokenMap, predecessor, disc, finish, time);
            TVertexIterator it_label(g);
            for(;!atEnd(it_label);goNext(it_label)) {
                if ((getProperty(tokenMap, getValue(it_label)) == true) &&
                    (getProperty(components, getValue(it_label)) == nilPred)) {
                    assignProperty(components, getValue(it_label), label);
                }
            }
            ++label;
        }
    }

    return label;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_STRONGLY_CONNECTED_COMPNENTS_H_
