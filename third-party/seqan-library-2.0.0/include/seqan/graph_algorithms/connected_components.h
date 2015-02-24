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
// Implementation of the Connected-Components algorithm.
// ==========================================================================

#ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_CONNECTED_COMPONENTS_H_
#define INCLUDE_SEQAN_GRAPH_ALGORITHMS_CONNECTED_COMPONENTS_H_

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
// Function connectedComponents()
// ----------------------------------------------------------------------------

/*!
 * @fn connectedComponents
 *
 * @headerfile <seqan/graph_algorithms.h>
 *
 * @brief Decomposes an undirected graph into its connected components.
 *
 * @signature TSize connectedComponents(components, g);
 *
 * @param[out] components
 *               A property map.  Each vertex is mapped to a component id.  If two vertices share the same id
 *               they are in the same component.
 * @param[in]  g An undirected graph. Types: Undirected Graph
 *
 * @return TSize The number of components.
 *
 * @section Examples
 *
 * A simple example on how to use this function.
 *
 * @code{.cpp}
 * // Build Input.
 * Graph<Undirected<> > graph;
 * for (unsigned i = 0; i < 5; ++i)
 *     addVertex(graph);
 * addEdge(graph, 0, 1);
 * addEdge(graph, 0, 3);
 * addEdge(graph, 2, 4);
 * String<unsigned> components;
 * unsigned numComponents = 0;
 *
 * // Call Algorithm.
 * numComponents = connectedComponents(g, components);
 *
 * // Print Result.
 * std::cout << "Number of components: " << numComponents << std::endl;
 * std::cout << std::endl << "Vertex -> Component" << std::endl;
 * for (unsigned i = 0; i < length(components); ++i)
 *     std::cout << i << " -> " << components[i] << std::endl;
 * @endcode
 * The output now is:
 *
 * @code{.console}
 * Number of components: 2
 *
 * Vertex -> Component
 * 0 -> 0
 * 1 -> 0
 * 2 -> 1
 * 3 -> 0
 * 4 -> 1
 * @endcode
 */
template <typename TSpec, typename TComponents>
typename Size<Graph<TSpec> >::Type
connectedComponents(TComponents & components,
                    Graph<TSpec> const & g)
{
    typedef typename VertexDescriptor<Graph<TSpec> >::Type TVertexDescriptor;

    // Initialize Union-Find data structure.
    UnionFind<TVertexDescriptor> uf;
    resizeVertexMap(uf, g);

    // Use Union-Find data structure to compute connected components.
    for (typename Iterator<Graph<TSpec>, EdgeIterator>::Type it(g); !atEnd(it); goNext(it))
        joinSets(uf, findSet(uf, sourceVertex(it)), findSet(uf, targetVertex(it)));

    // Build final component map.
    resizeVertexMap(components, g);
    unsigned c = 0;
    std::map<TVertexDescriptor, unsigned> reprToComponent;
    for (typename Iterator<Graph<TSpec>, VertexIterator>::Type it(g); !atEnd(it); goNext(it))
    {
        TVertexDescriptor v = findSet(uf, *it);

        // Build mapping from representant to component id on the fly.
        if (!reprToComponent.count(v))
            reprToComponent[v] = c++;

        assignProperty(components, *it, reprToComponent[v]);
    }

    return c;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_CONNECTED_COMPONENTS_H_
