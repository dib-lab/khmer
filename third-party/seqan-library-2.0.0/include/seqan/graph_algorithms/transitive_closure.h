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
// Implementation of Transitive-Closure algorithm.
// ==========================================================================

#ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_TRANSITIVE_CLOSURE_H_
#define INCLUDE_SEQAN_GRAPH_ALGORITHMS_TRANSITIVE_CLOSURE_H_

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
// Function transitiveClosure()
// ----------------------------------------------------------------------------

/*!
 * @fn transitiveClosure
 * @headerfile <seqan/graph_algorithms.h>
 * @brief Determines whether there is a path between any two given vertices or not.
 *
 * @signature void transitiveClosure(closure, g);
 *
 * @param[out] closure A matrix which indicates the closure.  Entry (i,j) in this matrix indicates whether there is a
 *                     path from i to j in the graph or not.  Types: Matrix
 * @param[in]  g       A directed graph.  Types: Directed Graph
 *
 * @section Example
 *
 * @include demos/graph_algorithms/transitive_closure.cpp
 *
 * @include demos/graph_algorithms/transitive_closure.cpp.stdout
 */
template <typename TSpec, typename TMatrix>
void transitiveClosure(TMatrix & closure,
                       Graph<TSpec> const & g)
{
    typedef typename Size<TMatrix>::Type TSize;

    // Initialize first closure matrix
    getAdjacencyMatrix(g,closure);
    TSize len = (TSize) std::sqrt((double) length(closure));
    for (TSize diag=0;diag < len;++diag) assignValue(closure, diag*len+diag,1);

    // Transitive Closure
    TMatrix local = closure;
    for (TSize k=0;k<len;++k) {
        for(TSize i=0;i<len;++i) {
            for(TSize j=0;j<len;++j) {
                bool t_ij = static_cast<int>(getValue(closure, i*len+j)) > 0;
                bool t_ik = static_cast<int>(getValue(closure, i*len+k)) > 0;
                bool t_kj = static_cast<int>(getValue(closure, k*len+j)) > 0;
                assignValue(local, i*len+j, t_ij || (t_ik && t_kj));
            }
        }
        closure = local;
    }
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_TRANSITIVE_CLOSURE_H_
