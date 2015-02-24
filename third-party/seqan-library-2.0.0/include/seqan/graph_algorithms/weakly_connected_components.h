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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Implementation of Weakly-Connected-Components algorithm.
// ==========================================================================

#ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_WEAKLY_CONNECTED_COMPONENTS_H_
#define INCLUDE_SEQAN_GRAPH_ALGORITHMS_WEAKLY_CONNECTED_COMPONENTS_H_

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
// Function weaklyConnectedComponents()
// ----------------------------------------------------------------------------

/*!
 * @fn weaklyConnectedComponents
 * @headerfile <seqan/graph_algorithms.h>
 * @brief Compute weakly connected components of a directed graph.
 *
 * @signature TSize weaklyConnectedComponents(components, g);
 *
 * @param[out] components
 *               A property map.  Each vertex is mapped to a component id.  If two vertices share the same id they
 *               are in the same component.
 * @param[in]  g A @link DirectedGraph @endlink to use for the input.
 *
 * @return TSize The number of weakly connected components  (Metafunction: @link Graph#Size @endlink of the type
 *               of <tt>g</tt>).
 *
 * The running time is <tt>O(n a(n, n))</tt> where <tt>a</tt> is the inverse Ackermann function and thus almost linear.
 * The union find data structure is used since the graph implementations do not allow the efficient iteration of
 * in-edges.
 */
template <typename TSpec, typename TComponents>
typename Size<Graph<TSpec> >::Type
weaklyConnectedComponents(TComponents & components,
                          Graph<TSpec> const & g)

{
    typedef Graph<TSpec> TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

    // Initialization.
    UnionFind<TVertexDescriptor> unionFind;
    resizeVertexMap(unionFind, g);

    // Iterate over all edges, joining weakly connected components.
    for (TEdgeIterator itE(g); !atEnd(itE); goNext(itE))
        joinSets(unionFind, findSet(unionFind, sourceVertex(itE)), findSet(unionFind, targetVertex(itE)));

    // Count number of sets.
    TSize setCount = 0;
    for (TVertexIterator itV(g); !atEnd(itV); goNext(itV))
        setCount += (findSet(unionFind, *itV) == *itV);

    // Build a map from graph vertex descriptor to component id.
    TSize nextId = 0;
    clear(components);
    resizeVertexMap(components, g, setCount);  // setCount is sentinel value
    for (TVertexIterator itV(g); !atEnd(itV); goNext(itV)) {
        if (getProperty(components, findSet(unionFind, *itV)) == setCount)
            assignProperty(components, findSet(unionFind, *itV), nextId++);
        assignProperty(components, *itV, getProperty(components, findSet(unionFind, *itV)));
    }

    return setCount;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_WEAKLY_CONNECTED_COMPONENTS_H_
