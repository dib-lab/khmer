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
// Implementation of Topological-Sort algorithm.
// ==========================================================================

#ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_TOPOLOGICAL_SORT_H_
#define INCLUDE_SEQAN_GRAPH_ALGORITHMS_TOPOLOGICAL_SORT_H_

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
// Function topologicalSort()
// ----------------------------------------------------------------------------

/*!
 * @fn topologicalSort
 * @headerfile <seqan/graph_algorithms.h>
 * @brief Performs a topological sort on a directed acyclic graph (DAG).
 *
 * @signature void topologicalSort(topSort, g);
 *
 * @param[out] topSort A topological ordering of the vertices. Types: String.
 * @param[in]  g       A directed acyclic graph. Types: Directed Graph
 *
 * A topological sort is a linear ordering of all its vertices such that if the graph contains an edge (u,v) then u
 * appears before v in the ordering.
 *
 * @section Example
 *
 * @include demos/graph_algorithms/topological_sort.cpp
 *
 * @include demos/graph_algorithms/topological_sort.cpp.stdout
 */
template <typename TSpec, typename TVertexDescriptor>
void topologicalSort(String<TVertexDescriptor> & topSort,
                     Graph<TSpec> const & g)
{
    typedef typename Size<Graph<TSpec> >::Type TSize;

    // Variable definition.
    String<TSize> predMap;
    String<TSize> discoveryTimeMap;
    String<TSize> finishingTimeMap;

    // Perform DFS.
    depthFirstSearch(predMap, discoveryTimeMap, finishingTimeMap, g);
    SEQAN_ASSERT_EQ(numVertices(g), length(predMap));
    SEQAN_ASSERT_EQ(numVertices(g), length(discoveryTimeMap));
    SEQAN_ASSERT_EQ(numVertices(g), length(finishingTimeMap));

    // Order vertices.
    typedef std::pair<TSize, TVertexDescriptor> TTimeVertexPair;
    std::priority_queue<TTimeVertexPair> q;
    typedef typename Iterator<Graph<TSpec>, VertexIterator>::Type TVertexIterator;
    TVertexIterator it(g);
    for (; !atEnd(it); goNext(it))
        q.push(std::make_pair(getProperty(finishingTimeMap, getValue(it)), getValue(it)));

    // Create topological order.
    resize(topSort, numVertices(g));
    TSize count = 0;
    while (!q.empty())
    {
        assignValue(topSort, count, q.top().second);
        q.pop();
        ++count;
    }
    SEQAN_ASSERT_EQ(length(topSort), numVertices(g));
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_TOPOLOGICAL_SORT_H_
