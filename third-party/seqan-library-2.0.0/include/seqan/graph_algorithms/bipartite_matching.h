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
// Implementation of the Bipartite Matching algorithm.
//
// WARNING: Functionality not carefully tested!
// ==========================================================================

#ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_BIPARTITE_MATCHING_H_
#define INCLUDE_SEQAN_GRAPH_ALGORITHMS_BIPARTITE_MATCHING_H_

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
// Function bipartiteMatching()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TVertexMap, typename TEdges>
typename Size<Graph<TSpec> >::Type
bipartiteMatching(String<TEdges> & edges,
                  Graph<TSpec> const & g,
                  TVertexMap const & vertMap)
{
    typedef Graph<TSpec> TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIter;

    clear(edges);
    TVertexDescriptor source = addVertex(g);
    TVertexDescriptor target = addVertex(g);
    TVertexIter itV(g);
    for(;!atEnd(itV); goNext(itV)) {
        if ((value(itV) != source) && (value(itV) != target)) {
            if (getProperty(vertMap, value(itV)) == false) {
                addEdge(g, source, value(itV));
            } else {
                addEdge(g, value(itV), target);
            }
        }
    }

    // Use Ford-Fulkerson to determine a matching
    String<TSize> capMap;
    resizeEdgeMap(capMap, g);
    typedef typename Iterator<String<TSize> >::Type TCapIter;
    TCapIter capIt = begin(capMap);
    TCapIter capItEnd = end(capMap);
    for(;capIt != capItEnd; ++capIt) value(capIt) = 1;
    String<TSize> flow;
    TSize valF = fordFulkersonAlgorithm(g, source, target, capMap, flow);

    typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
    TEdgeIterator itEdge(g);
    for(;!atEnd(itEdge);goNext(itEdge)) {
        if (getProperty(flow, getValue(itEdge)) == 1) {
            TVertexDescriptor sV = sourceVertex(itEdge);
            TVertexDescriptor tV = targetVertex(itEdge);
            if ((sV != source) && (tV != target)) appendValue(edges, TEdges(sV, tV));
        }
    }
    removeVertex(g, source);
    removeVertex(g, target);

    return valF;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_BIPARTITE_MATCHING_H_
