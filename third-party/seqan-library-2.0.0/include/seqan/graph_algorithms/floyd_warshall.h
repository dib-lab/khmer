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
// Implementation of the Floyd-Warshall algorithm.
// ==========================================================================

#ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_FLOYD_WARSHALL_H_
#define INCLUDE_SEQAN_GRAPH_ALGORITHMS_FLOYD_WARSHALL_H_

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
// Function floydWarshallAlgorithm()
// ----------------------------------------------------------------------------

/*!
 * @fn floydWarshallAlgorithm
 * @headerfile <seqan/graph_algorithms.h>
 * @brief Finds shortest paths between all pairs of vertices in a graph.
 *
 * @signature void floydWarshallAlgorithm(g, weight, distance, predecessor);
 *
 * @param[out] predecessor A matrix with predecessors.  Entry (i,j) in this matrix indicates the predecessor of j on
 *                         a shortest path from vertex i to vertex j.  You can use <tt>_printAllPairsShortestPath(g,
 *                         predecessor, i, j)</tt> to print the shortest path from i to j.  Types: Matrix
 * @param[out] distance    A matrix with distances.Entry (i,j) in this matrix indicates the distance from vertex i
 *                         to vertex j.  Types: Matrix
 * @param[in]  weight      A weight map.  A property map with edge weights.  Edge weights may be negative.
 * @param[in]  g           A directed graph.  Types: Directed Graph
 *
 * The graph must be free of negative-weight cycles.
 *
 * @section Example
 *
 * @include demos/graph_algorithms/floyd_warshall_algorithm.cpp
 *
 * @include demos/graph_algorithms/floyd_warshall_algorithm.cpp.stdout
 *
 * @see allPairsShortestPath
 */
template <typename TSpec, typename TWeightMap, typename TMatrix, typename TPredecessor>
void floydWarshallAlgorithm(TMatrix & distMatrix,
                            TPredecessor & predecessor,
                            Graph<TSpec> const & g,
                            TWeightMap const & weight)
{
    typedef typename Size<TMatrix>::Type TSize;
    typedef typename Value<TMatrix>::Type TMatrixVal;

    // Initialize first distance matrix
    _initializeAllPairs(g,weight,distMatrix,predecessor);

    // Floyd-Warshall
    TSize len = (TSize) std::sqrt((double) length(distMatrix));
    TMatrix local = distMatrix;
    for(TSize k=0;k<len;++k) {
        for(TSize i=0;i<len;++i) {
            for(TSize j=0;j<len;++j) {
                TMatrixVal min1 = getValue(distMatrix, i*len+j);
                TMatrixVal min2 = getValue(distMatrix, i*len+k) + getValue(distMatrix, k*len + j);
                if (min2 < min1) {
                    assignValue(local, i*len+j,min2);
                    assignValue(predecessor, i*len+j,getValue(predecessor, k*len+j));
                } else {
                    assignValue(local, i*len+j,min1);
                    assignValue(predecessor, i*len+j, getValue(predecessor, i*len+j));
                }
            }
        }
        distMatrix=local;
    }
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_FLOYD_WARSHALL_H_
