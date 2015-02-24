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
// Implementation of All-Pairs-Shortest-Path algorithm.
// ==========================================================================

#ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_ALL_PAIRS_SHORTEST_PATH_H_
#define INCLUDE_SEQAN_GRAPH_ALGORITHMS_ALL_PAIRS_SHORTEST_PATH_H_

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
// Function allPairsShortestPath()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TPredecessor, typename TVertexDescriptor>
inline void
_printAllPairsShortestPath(Graph<TSpec> const& g,
                           TPredecessor& predecessor,
                           TVertexDescriptor const i,
                           TVertexDescriptor const j)
{
    typedef typename Size<TPredecessor>::Type TSize;
    TSize len = getIdUpperBound(g.data_id_managerV);
    if (i==j) {
        std::cout << i;
    } else if (getValue(predecessor, i*len+j) == getNil<typename VertexDescriptor<Graph<TSpec> >::Type>()) {
        std::cout << "No path from " << i << " to " << j << " exists.";
    } else {
        _printAllPairsShortestPath(g,predecessor, i, (TVertexDescriptor) getValue(predecessor, i*len+j));
        std::cout << "," << j;
    }
}

template <typename TSpec, typename TWeightMap, typename TMatrix, typename TPredecessor>
void
_initializeAllPairs(Graph<TSpec> const& g,
                        TWeightMap const& weight,
                        TMatrix& matrix,
                        TPredecessor& predecessor)
{
    typedef Graph<TSpec> TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
    typedef typename Size<TMatrix>::Type TSize;
    typedef typename Value<TWeightMap>::Type TWeightVal;
    typedef typename Value<TPredecessor>::Type TPredVal;

    // Create adjacency-like matrix
    TSize len = getIdUpperBound(g.data_id_managerV);
    resize(matrix, len * len);
    resize(predecessor, len * len);
    TWeightVal infWeight = _getInfinityDistance(weight);
    TPredVal nilPred = getNil<TVertexDescriptor>();
    for (TSize row=0;row < len;++row) {
        for (TSize col=0;col < len;++col) {
            if (row != col) assignValue(matrix, row*len + col, infWeight);
            else assignValue(matrix, row*len + col, 0);
            assignValue(predecessor, row*len + col, nilPred);
        }
    }

    // Include edge weights and initial predecessors
    TVertexIterator it(g);
    for(;!atEnd(it);goNext(it)) {
        TVertexDescriptor u = getValue(it);
        TOutEdgeIterator itout(g, u);
        for(;!atEnd(itout);++itout) {
            TVertexDescriptor v = targetVertex(g,getValue(itout));
            assignValue(matrix, u*len + v, getProperty(weight, getValue(itout)));
            assignValue(predecessor, u*len + v, u);
        }
    }
}

template <typename TMatrix, typename TPredecessor, typename TInfDist>
void
_extendShortestPaths(TMatrix& local,
                       TMatrix& w,
                       TPredecessor& predecessor,
                       TInfDist const infDist)
{
    typedef typename Value<TMatrix>::Type TMatrixVal;
    typedef typename Value<TPredecessor>::Type TPredVal;
    typedef typename Size<TMatrix>::Type TSize;
    TMatrix oldLocal = local;
    TPredecessor oldPredecessor = predecessor;
    TSize len = (TSize) std::sqrt((double) length(oldLocal));
    for(TSize i = 0; i<len;++i) {
        for(TSize j = 0; j<len;++j) {
            if (i==j) continue;
            assignValue(local, i*len+j,infDist);
            TPredVal ind = 0;
            for(TSize k = 0; k<len;++k) {
                TMatrixVal min1 = getValue(local, i*len+j);
                TMatrixVal min2 = getValue(oldLocal, i*len+k) + getValue(w, k*len + j);
                if (min2 < min1) {
                    assignValue(local, i*len+j,min2);
                    ind = k;
                }
            }
            if (getValue(oldLocal, i*len+j) > getValue(local, i*len+j)) {
                assignValue(predecessor, i*len+j,ind);
            }
        }
    }
}

/*!
 * @fn allPairsShortestPath
 * @headerfile <seqan/graph_algorithms.h>
 * @brief Finds shortest paths between all pairs of vertices in a graph.
 * @signature void allPairsShortestPath(graph, weight, distance, predecessor);
 *
 * @param[out] distance    A @link Matrix @endlink with distances. Entry <tt>(i,j)</tt> in this matrix indicates the
 *                         distance from vertex <tt>i</tt> to vertex <tt>j</tt>.
 * @param[out] predecessor A @link Matrix @endlink with predecessors. Entry <tt>(i,j)</tt> in this matrix indicates the
 *                         predecessor of <tt>j</tt> on a shortest path from vertex <tt>i</tt> to vertex <tt>j</tt>.
 *                         You can use <tt>_printAllPairsShortestPath(graph, predecessor, i, j)</tt> to print the
 *                         shortest path from <tt>i</tt> to <tt>j</tt>.
 * @param[in]  graph       A @link DirectedGraph Directed Graph @endlink.
 * @param[in]  weight      A property map with edge weights. Edge weights may be negative.
 *
 * @section Example
 *
 * @include demos/graph_algorithms/all_pairs_shortest_path.cpp
 *
 * @include demos/graph_algorithms/all_pairs_shortest_path.cpp.stdout
 *
 * @see floydWarshallAlgorithm
 */
template <typename TSpec, typename TWeightMap, typename TMatrix, typename TPredecessor>
void allPairsShortestPath(TMatrix & distMatrix,
                          TPredecessor & predecessor,
                          Graph<TSpec> const & g,
                          TWeightMap const & weight)
{
    typedef typename Size<TMatrix>::Type TSize;
    typedef typename Value<TWeightMap>::Type TWeightVal;
    TWeightVal infWeight = _getInfinityDistance(weight);

    // Initialize first distance matrix
    _initializeAllPairs(g,weight,distMatrix,predecessor);

    TSize len = (TSize) sqrt((double) length(distMatrix));
    TMatrix local = distMatrix;
    for(TSize m=2;m<len;++m) {
        _extendShortestPaths(local,distMatrix,predecessor, infWeight);
    }
    distMatrix = local;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_ALL_PAIRS_SHORTEST_PATH_H_
