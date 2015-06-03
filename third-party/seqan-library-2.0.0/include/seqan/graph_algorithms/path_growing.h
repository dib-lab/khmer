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
// Implementation of Path-Growing algorithm.
//
// WARNING: Functionality not carefully tested!
// ==========================================================================

#ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_PATH_GROWING_H_
#define INCLUDE_SEQAN_GRAPH_ALGORITHMS_PATH_GROWING_H_

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
// Functions
// ----------------------------------------------------------------------------

template <typename TSpec, typename TWeightMap, typename TEdgeMap>
typename Value<TWeightMap>::Type
pathGrowingAlgorithm(TEdgeMap & edgeMap1,
                     Graph<TSpec> const & g,
                     TWeightMap const & weightMap)
{
    typedef Graph<TSpec> TGraph;
    typedef typename Value<TWeightMap>::Type TValue;
    typedef typename Size<Graph<TSpec> >::Type TSize;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

    // Make a copy of the graph
    TGraph mutant(g);

    // Initialy not a single edge is selected
    resize(edgeMap1, getIdUpperBound(_getEdgeIdManager(g)), false);
    TEdgeMap edgeMap2 = edgeMap1;
    TValue edgeMap1Sum = 0;
    TValue edgeMap2Sum = 0;

    // Run the algorithm
    TSize i = 1;
    while (numEdges(mutant) > 0) {
        TVertexIterator itVert(mutant);
        while (outDegree(mutant, *itVert) < 1) goNext(itVert);
        TVertexDescriptor x = *itVert;
        TVertexDescriptor y;
        while (outDegree(mutant, x) >= 1) {
            TOutEdgeIterator itOut(mutant, x);
            TEdgeDescriptor e = *itOut;
            TValue max = getProperty(weightMap, e);
            y = targetVertex(itOut);
            goNext(itOut);
            for(;!atEnd(itOut);++itOut) {
                if (getProperty(weightMap, *itOut) > max) {
                    e = *itOut;
                    max = getProperty(weightMap, e);
                    y = targetVertex(itOut);
                }
            }
            if (i == 1) {
                // Mark the edge for m1
                assignProperty(edgeMap1, e, true);
                edgeMap1Sum += max;
            } else {
                // Mark the edge for m2
                assignProperty(edgeMap2, e, true);
                edgeMap2Sum += max;
            }
            i = 3 - i;
            removeVertex(mutant, x);
            x = y;
        }
    }


    // Check whether we have to swap bool arrays
    if (edgeMap2Sum > edgeMap1Sum) {
        edgeMap1Sum = edgeMap2Sum;
        edgeMap1 = edgeMap2;
    }

    return edgeMap1Sum;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_PATH_GROWING_H_
