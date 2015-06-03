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

#ifndef SEQAN_HEADER_GRAPH_BASE_H
#define SEQAN_HEADER_GRAPH_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// General Graph Metafunction
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/*!
 * @mfn Graph#EdgeDescriptor
 * @brief Type of an object that represents an edge descriptor.
 *
 * @signature EdgeDescriptor<T>::Type;
 *
 * @tparam T The @link Graph @endlink type to query.
 *
 * @return Type The resulting edge descriptor type.
 *
 * The edge descriptor is a unique handle to a given edge in a graph.  It is used in various graph functions, e.g. to
 * remove edges, to assign cargo to an edge or to get the end points of an edge.  It is also used to attache properties
 * to edges via property maps.
 *
 * @section Examples
 *
 * @code{.cpp}
 * EdgeDescriptor<Graph<> >::Type ed;  // eD is an edge descriptor
 * @endcode
 */

template<typename T>
struct EdgeDescriptor;

//////////////////////////////////////////////////////////////////////////////

/*!
 * @mfn Graph#Cargo
 * @brief Type for the graph's cargo.
 *
 * @signature Cargo<T>::Type;
 *
 * @tparam T The @link Graph @endlink type to query.
 *
 * @return Type The resulting cargo type.
 */

template<typename T>
struct Cargo;


//////////////////////////////////////////////////////////////////////////////

/*!
 * @mfn Graph#EdgeType
 * @brief Edge type of a graph class.
 *
 * @signature EdgeType<T>::Type;
 *
 * @tparam T The @link Graph @endlink type to query.
 *
 * @return Type The resulting edge stump type that is used in the graph.
 *
 * @section Examples
 *
 * @code{.cpp}
 * EdgeType<TGraph>::Type e;  // e is an edge in a TGraph
 * @endcode
 */

template<typename T>
struct EdgeType;

//////////////////////////////////////////////////////////////////////////////

/*!
 * @concept GraphOverAlphabetConcept
 * @brief A graph construted over an alphabet.
 *
 * @mfn GraphOverAlphabetConcept#Alphabet
 * @brief Return the Alphabe type of a graph over an alphabet.
 *
 * @signature Alphabet<T>::Type;
 *
 * @tparam T The @link GraphOverAlphabetConcept @endlink type to query.
 *
 * @return Type The alphabe type.
 *
 * @section Examples
 *
 * @code{.cpp}
 * Alphabet<Graph<Automaton<Dna> > >::Type c;  // c is of type Dna
 * @endcode
 */

template<typename T>
struct Alphabet;


//////////////////////////////////////////////////////////////////////////////

/*!
 * @mfn Graph#EdgeIdHandler
 * @brief Type of an object that represents an IdManager.
 *
 * @signature EdgeIdHandler<T>::Type;
 *
 * @tparam T The Graph to query.
 *
 * @return Type The IdManager type.
 *
 * The exact IdManager type depends on the edge stump type.  If the edge stump has no ids then the IdManager simply
 * counts edge ids and otherwise it manages a list of free and used ids.
 */

template<typename T>
struct EdgeIdHandler;


//////////////////////////////////////////////////////////////////////////////

/*!
 * @mfn Graph#VertexIdHandler
 * @brief Type of an object that repreestns an IdManager.
 *
 * @signature VertexIdHandler<T>::Type;
 *
 * @tparam T The Graph to query.
 *
 * @return Type The IdManager type.
 */

template<typename T>
struct VertexIdHandler;


//////////////////////////////////////////////////////////////////////////////
// General Graph Tags
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

struct WithoutEdgeId_;
typedef Tag<WithoutEdgeId_> const WithoutEdgeId;

//////////////////////////////////////////////////////////////////////////////

struct TreeTag_;
typedef Tag<TreeTag_> const TreeTag;


//////////////////////////////////////////////////////////////////////////////
// Graph Iterator Tags
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/*!
 * @defgroup GraphIteratorTags Graph Iterator Tags
 * @brief Tags that can be used to get iterators on graphs.
 *
 * @tag GraphIteratorTags#VertexIterator
 * @brief Iterate over all vertices of a graph.
 *
 * @signature typedef Tag<VertexIterator_> const VertexIterator;
 *
 * @tag GraphIteratorTags#EdgeIterator
 * @brief Iterate over all edges of a graph.
 *
 * @signature typedef Tag<EdgeIterator_> const EdgeIterator;
 *
 * @tag GraphIteratorTags#OutEdgeIterator
 * @brief Iterate over all out edges of a vertex.
 *
 * @signature typedef Tag<OutEdgeIterator_> const OutEdgeIterator;
 *
 * @tag GraphIteratorTags#AdjacencyIterator
 * @brief Iterate over all adjacent vertices of a given vertex.
 *
 * @signature typedef Tag<AdjacencyIterator_> const AdjacencyIterator;
 *
 * @tag GraphIteratorTags#BfsIterator
 * @brief Iterate over all vertices of a graph in breadth-first fashion starting from a given vertex.
 *
 * @signature typedef Tag<BfsIterator_> const BfsIterator;
 *
 * @tag GraphIteratorTags#DfsPreorder
 * @brief Iterate over all vertices of a graph in depth-first fashion.
 *
 * @signature typedef Tag<DfsPreorder_> const DfsPreorder;
 */

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
struct VertexIterator_;
typedef Tag<VertexIterator_> const VertexIterator;

//////////////////////////////////////////////////////////////////////////////
struct EdgeIterator_;
typedef Tag<EdgeIterator_> const EdgeIterator;

//////////////////////////////////////////////////////////////////////////////
struct OutEdgeIterator_;
typedef Tag<OutEdgeIterator_> const OutEdgeIterator;

//////////////////////////////////////////////////////////////////////////////
struct AdjacencyIterator_;
typedef Tag<AdjacencyIterator_> const AdjacencyIterator;

//////////////////////////////////////////////////////////////////////////////
struct BfsIterator_;
typedef Tag<BfsIterator_> const BfsIterator;

//////////////////////////////////////////////////////////////////////////////
struct DfsPreorder_;
typedef Tag<DfsPreorder_> const DfsPreorder;





//////////////////////////////////////////////////////////////////////////////
// Graph - Default edge stump
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo = void, bool TList = true, bool TSource = false, bool TId = true, typename TSpec = Default>
class EdgeStump;

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
struct VertexDescriptor<EdgeStump<TCargo, TList, TSource, TId, TSpec> >
{
    typedef typename Id<EdgeStump<TCargo, TList, TSource, TId, TSpec> >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
struct VertexDescriptor<EdgeStump<TCargo, TList, TSource, TId, TSpec> const>
{
    typedef typename Id<EdgeStump<TCargo, TList, TSource, TId, TSpec> >::Type Type;
};




//////////////////////////////////////////////////////////////////////////////
// Graph - Default Id Manager
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TIdType = unsigned int, typename TSpec = Default>
class IdManager;

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, typename TSpec>
struct EdgeIdHandler<EdgeStump<TCargo, TList, TSource, false, TSpec> > {
    typedef IdManager<void> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, typename TSpec>
struct EdgeIdHandler<EdgeStump<TCargo, TList, TSource, true, TSpec> > {
    typedef IdManager<typename Id<EdgeStump<TCargo, TList, TSource, true, TSpec> >::Type> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename T>
struct VertexIdHandler {
    typedef IdManager<> Type;
};

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
