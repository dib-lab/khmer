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

#ifndef SEQAN_HEADER_GRAPH_EDGESTUMP_H
#define SEQAN_HEADER_GRAPH_EDGESTUMP_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
//    Graph - EdgeStump
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/*!
 * @class EdgeStump
 * @headerfile <seqan/graph_types.h>
 * @brief Encapsulate information for a single edge.
 *
 * The EdgeStump either represents a list entry in the adjacency of a graph or an array field if edges are stored in an
 * array.
 *
 * @signature template <typename TCargo, bool IS_LIST, bool STORE_SOURCE, bool STORE_ID, typename TSpec>
 *            class EdgeStump;
 *
 * @tparam TCargo       The cargo type of an edge.  The cargo can be used to store arbitrary information with an
 *                      edge or be <tt>void</tt>.  Default: <tt>void</tt>.
 * @tparam IS_LIST      A bool value that indicates whether it is a list or not.  Default: <tt>true</tt>.
 * @tparam STORE_SOURCE A bool value that indicates whether the source is stored in the EdgeStump or not.
 *                      Default: <tt>false</tt>.
 * @tparam STORE_ID     A bool value that indictes whether the id is tored in the EdgeStump or not.  Note: Without edge
 *                      ids, external property maps do not work for edges!  Default: <tt>true</tt>.
 * @tparam TSpec        The specializing type.  Default: <tt>Default</tt>.
 *
 * @section Remarks
 *
 * The default EdgeStump in all graph types does not consider a cargo.  However, ni default usage every graph does store
 * an edge id.  Edge ids are used to append additional properties to edges with the help of external property maps.
 */

template<typename TCargo, typename TSpec>
class EdgeStump<TCargo, true, false, false, TSpec>
{
public:
    EdgeStump() : data_target(), data_cargo(), data_nextT()
    {}

    typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor_;
    TVertexDescriptor_ data_target;
    TCargo data_cargo;
    EdgeStump* data_nextT;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
class EdgeStump<TCargo, true, false, true, TSpec>
{
public:
    EdgeStump() : data_target(), data_id(), data_cargo(), data_nextT()
    {}

    typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor_;
    typedef typename Id<EdgeStump>::Type TId_;
    TVertexDescriptor_ data_target;
    TId_ data_id;
    TCargo data_cargo;
    EdgeStump* data_nextT;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
class EdgeStump<TCargo, true, true, false, TSpec>
{
public:
    EdgeStump() : data_target(), data_source(), data_cargo(), data_nextT(), data_nextS()
    {}

    typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor_;
    TVertexDescriptor_ data_target;
    TVertexDescriptor_ data_source;
    TCargo data_cargo;
    EdgeStump* data_nextT;
    EdgeStump* data_nextS;
};

//////////////////////////////////////////////////////////////////////////////


template<typename TCargo, typename TSpec>
class EdgeStump<TCargo, true, true, true, TSpec>
{
public:
    EdgeStump() : data_target(), data_source(), data_id(), data_cargo(), data_nextT(), data_nextS()
    {}

    typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor_;
    typedef typename Id<EdgeStump>::Type TId_;
    TVertexDescriptor_ data_target;
    TVertexDescriptor_ data_source;
    TId_ data_id;
    TCargo data_cargo;
    EdgeStump* data_nextT;
    EdgeStump* data_nextS;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
class EdgeStump<TCargo, false, false, false, TSpec>
{
public:
    EdgeStump() : data_target(), data_cargo()
    {}

    typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor_;
    TVertexDescriptor_ data_target;
    TCargo data_cargo;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
class EdgeStump<TCargo, false, false, true, TSpec>
{
public:
    EdgeStump() : data_target(), data_id(), data_cargo()
    {}

    typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor_;
    typedef typename Id<EdgeStump>::Type TId_;
    TVertexDescriptor_ data_target;
    TId_ data_id;
    TCargo data_cargo;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
class EdgeStump<TCargo, false, true, false, TSpec>
{
public:
    EdgeStump() : data_target(), data_source(), data_cargo()
    {}

    typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor_;
    TVertexDescriptor_ data_target;
    TVertexDescriptor_ data_source;
    TCargo data_cargo;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
class EdgeStump<TCargo, false, true, true, TSpec>
{
public:
    EdgeStump() : data_target(), data_source(), data_id(), data_cargo()
    {}

    typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor_;
    typedef typename Id<EdgeStump>::Type TId_;
    TVertexDescriptor_ data_target;
    TVertexDescriptor_ data_source;
    TId_ data_id;
    TCargo data_cargo;
};

//////////////////////////////////////////////////////////////////////////////
//    Graph - Cargoless EdgeStump
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
class EdgeStump<void, true, false, false, TSpec>
{
public:
    EdgeStump() : data_target(), data_nextT()
    {}

    typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor_;
    TVertexDescriptor_ data_target;
    EdgeStump* data_nextT;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
class EdgeStump<void, true, false, true, TSpec>
{
public:
    EdgeStump() : data_target(), data_id(), data_nextT()
    {}

    typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor_;
    typedef typename Id<EdgeStump>::Type TId_;
    TVertexDescriptor_ data_target;
    TId_ data_id;
    EdgeStump* data_nextT;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
class EdgeStump<void, true, true, false, TSpec>
{
public:
    EdgeStump() : data_target(), data_source(), data_nextT(), data_nextS()
    {}

    typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor_;
    TVertexDescriptor_ data_target;
    TVertexDescriptor_ data_source;
    EdgeStump* data_nextT;
    EdgeStump* data_nextS;
};

//////////////////////////////////////////////////////////////////////////////


template<typename TSpec>
class EdgeStump<void, true, true, true, TSpec>
{
public:
    EdgeStump() : data_target(), data_source(), data_id(), data_nextT(), data_nextS()
    {}

    typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor_;
    typedef typename Id<EdgeStump>::Type TId_;
    TVertexDescriptor_ data_target;
    TVertexDescriptor_ data_source;
    TId_ data_id;
    EdgeStump* data_nextT;
    EdgeStump* data_nextS;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
class EdgeStump<void, false, false, false, TSpec>
{
public:
    EdgeStump() : data_target()
    {}

    typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor_;
    TVertexDescriptor_ data_target;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
class EdgeStump<void, false, false, true, TSpec>
{
public:
    EdgeStump() : data_target(), data_id()
    {}

    typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor_;
    typedef typename Id<EdgeStump>::Type TId_;
    TVertexDescriptor_ data_target;
    TId_ data_id;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
class EdgeStump<void, false, true, false, TSpec>
{
public:
    EdgeStump() : data_target(), data_source()
    {}

    typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor_;
    TVertexDescriptor_ data_target;
    TVertexDescriptor_ data_source;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
class EdgeStump<void, false, true, true, TSpec>
{
public:
    EdgeStump() : data_target(), data_source(), data_id()
    {}

    typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor_;
    typedef typename Id<EdgeStump>::Type TId_;
    TVertexDescriptor_ data_target;
    TVertexDescriptor_ data_source;
    TId_ data_id;
};

//////////////////////////////////////////////////////////////////////////////
// EdgeStump - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
struct Cargo<EdgeStump<TCargo, TList, TSource, TId, TSpec> > {
    typedef TCargo Type;
};

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
struct Cargo<EdgeStump<TCargo, TList, TSource, TId, TSpec> const> {
    typedef TCargo const Type;
};


template<bool TList, bool TSource, bool TId, typename TSpec>
struct Cargo<EdgeStump<void, TList, TSource, TId, TSpec> > {
    typedef void* Type;
};

template<bool TList, bool TSource, bool TId, typename TSpec>
struct Cargo<EdgeStump<void, TList, TSource, TId, TSpec> const> {
    typedef void* Type;
};


//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
struct Spec<EdgeStump<TCargo, TList, TSource, TId, TSpec> >
{
    typedef TSpec Type;
};

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
struct Spec<EdgeStump<TCargo, TList, TSource, TId, TSpec> const>
{
    typedef TSpec Type;
};

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn EdgeStump#getCargo
 * @brief Return cargo for an EdgeStump.
 *
 * @signature TCargo getCargo(stump);
 *
 * @param[in] stump Pointer to the EdgeStump to query for its cargo.
 *
 * @return TCargo Reference to the cargo of the EdgeStump.
 */


template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
inline typename Cargo<EdgeStump<TCargo, TList, TSource, TId, TSpec> const>::Type&
getCargo(EdgeStump<TCargo, TList, TSource, TId, TSpec> const* es)
{
    SEQAN_CHECKPOINT
    return es->data_cargo;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
inline typename Cargo<EdgeStump<TCargo, TList, TSource, TId, TSpec> >::Type&
getCargo(EdgeStump<TCargo, TList, TSource, TId, TSpec>* es)
{
    SEQAN_CHECKPOINT
    return es->data_cargo;
}

//////////////////////////////////////////////////////////////////////////////

template<bool TList, bool TSource, bool TId, typename TSpec>
inline typename Cargo<EdgeStump<void, TList, TSource, TId, TSpec> const>::Type
getCargo(EdgeStump<void, TList, TSource, TId, TSpec> const*)
{
    SEQAN_CHECKPOINT
    // No real cargo
    return 0;
}

//////////////////////////////////////////////////////////////////////////////

template<bool TList, bool TSource, bool TId, typename TSpec>
inline typename Cargo<EdgeStump<void, TList, TSource, TId, TSpec> >::Type
getCargo(EdgeStump<void, TList, TSource, TId, TSpec>*)
{
    SEQAN_CHECKPOINT
    // No real cargo
    return 0;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn EdgeStump#cargo
 * @brief Return cargo for an EdgeStump.
 *
 * @signature TCargo cargo(stump);
 *
 * @param[in] stump Pointer to the EdgeStump to query for its cargo.
 *
 * @return TCargo Reference to the cargo of the EdgeStump.
 */


template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
inline typename Cargo<EdgeStump<TCargo, TList, TSource, TId, TSpec> const>::Type&
cargo(EdgeStump<TCargo, TList, TSource, TId, TSpec> const* es)
{
    SEQAN_CHECKPOINT
    return es->data_cargo;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
inline typename Cargo<EdgeStump<TCargo, TList, TSource, TId, TSpec> >::Type&
cargo(EdgeStump<TCargo, TList, TSource, TId, TSpec>* es)
{
    SEQAN_CHECKPOINT
    return es->data_cargo;
}


//////////////////////////////////////////////////////////////////////////////

template<bool TList, bool TSource, bool TId, typename TSpec>
inline typename Cargo<EdgeStump<void, TList, TSource, TId, TSpec> >::Type
cargo(EdgeStump<void, TList, TSource, TId, TSpec>*)
{
    SEQAN_CHECKPOINT
    // No real cargo
    return 0;
}

//////////////////////////////////////////////////////////////////////////////

template<bool TList, bool TSource, bool TId, typename TSpec>
inline typename Cargo<EdgeStump<void, TList, TSource, TId, TSpec> const>::Type
cargo(EdgeStump<void, TList, TSource, TId, TSpec> const*)
{
    SEQAN_CHECKPOINT
    // No real cargo
    return 0;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn EdgeStump#assignCargo
 * @brief Assigns a new cargo to the edge.
 *
 * @signature void assignCargo(stump, cargo);
 *
 * @param[in] stump Pointer to the EdgeStump to set the cargo of.
 *
 * @return TCargo Reference to the cargo of the EdgeStump.
 *
 * Calling assignCargo on EdgeStump objects without cargo does nothing.
 */

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec, typename TCargo2>
inline void
assignCargo(EdgeStump<TCargo, TList, TSource, TId, TSpec>* es,
            TCargo2 const& t)
{
    SEQAN_CHECKPOINT
    es->data_cargo =  (TCargo) t;
}

//////////////////////////////////////////////////////////////////////////////

template<bool TList, bool TSource, bool TId, typename TSpec, typename TCargo2>
inline void
assignCargo(EdgeStump<void, TList, TSource, TId, TSpec>*,
            TCargo2 const&)
{
    SEQAN_CHECKPOINT
    // No real cargo
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn EdgeStump#assignTarget
 * @brief Assigns a target vertex to an edge.
 *
 * @signature void assignTarget(stump, t);
 *
 * @param[in,out] stump Pointer to the EdgeStump.
 * @param[in]     t     Vertex descriptor to assign as the target.
 */


template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec, typename TVertexDescriptor>
inline void
assignTarget(EdgeStump<TCargo, TList, TSource, TId, TSpec>* es,
             TVertexDescriptor const t)
{
    SEQAN_CHECKPOINT
    es->data_target = t;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn EdgeStump#target
 * @brief Access to the target of an EdgeStump.
 *
 * @signature TVertexDescriptor target(stump);
 *
 * @param[in] stump Pointer to the EdgeStump to access the target of.
 *
 * @return TVertexDescriptor Reference to the target vertex descriptor of stump.
 */

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TList, TSource, TId, TSpec> >::Type&
target(EdgeStump<TCargo, TList, TSource, TId, TSpec>* es)
{
    SEQAN_CHECKPOINT
    return es->data_target;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TList, TSource, TId, TSpec> >::Type
target(EdgeStump<TCargo, TList, TSource, TId, TSpec> const* es)
{
    SEQAN_CHECKPOINT
    return es->data_target;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn EdgeStump#getTarget
 * @brief Get method for the target.
 *
 * @signature TVertexDescriptor getTarget(stump);
 *
 * @param[in] stump Pointer to the EdgeStump to get the target of.
 *
 * @return TVertexDescriptor The vetex descriptor stored in stump.
 */

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TList, TSource, TId, TSpec> const>::Type
getTarget(EdgeStump<TCargo, TList, TSource, TId, TSpec> const* es)
{
    SEQAN_CHECKPOINT
    return es->data_target;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TList, TSource, TId, TSpec> >::Type
getTarget(EdgeStump<TCargo, TList, TSource, TId, TSpec>* es)
{
    SEQAN_CHECKPOINT
    return es->data_target;
}


//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn EdgeStump#assignSource
 * @brief Assigns a source vertex to an edge.
 *
 * @signature void assignSource(stump, t);
 *
 * @param[in,out] stump Pointer to the EdgeStump.
 * @param[in]     t     Vertex descriptor to assign as the source.
 *
 * @section Remarks
 *
 * A source vertex is not required in an edge stump.  However, EdgeStump objects can be configured to contain a source
 * vertex, as in undirected graphs.
 */

template<typename TCargo, bool TList, bool TId, typename TSpec, typename TVertexDescriptor>
inline void
assignSource(EdgeStump<TCargo, TList, true, TId, TSpec>* es,
             TVertexDescriptor const s)
{
    SEQAN_CHECKPOINT
    es->data_source = s;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TId, typename TSpec, typename TVertexDescriptor>
inline void
assignSource(EdgeStump<TCargo, TList, false, TId, TSpec>*,
             TVertexDescriptor const)
{
    SEQAN_CHECKPOINT
    // NOP
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn EdgeStump#source
 * @brief Access to the source of an EdgeStump.
 *
 * @signature TVertexDescriptor source(stump);
 *
 * @param[in] stump Pointer to the EdgeStump to access the source of.
 *
 * @return TVertexDescriptor Reference to the source vertex descriptor of stump.
 *
 * @section Remarks
 *
 * A source vertex is not required in an edge stump.  However, EdgeStump objects can be configured to contain a source
 * vertex, as in undirected graphs.
 */

template<typename TCargo, bool TList, bool TId, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TList, true, TId, TSpec> >::Type&
source(EdgeStump<TCargo, TList, true, TId, TSpec>* es)
{
    SEQAN_CHECKPOINT
    return es->data_source;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TId, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TList, true, TId, TSpec> >::Type
source(EdgeStump<TCargo, TList, true, TId, TSpec> const* es)
{
    SEQAN_CHECKPOINT
    return es->data_source;
}

//////////////////////////////////////////////////////////////////////////////


template<typename TCargo, bool TList, bool TId, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TList, false, TId, TSpec> >::Type
source(EdgeStump<TCargo, TList, false, TId, TSpec>*)
{
    SEQAN_CHECKPOINT
    // No source available
    return 0;
}

//////////////////////////////////////////////////////////////////////////////


template<typename TCargo, bool TList, bool TId, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TList, false, TId, TSpec> >::Type
source(EdgeStump<TCargo, TList, false, TId, TSpec> const*)
{
    SEQAN_CHECKPOINT
    // No source available
    return 0;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn EdgeStump#getSource
 * @brief Get method for the source.
 *
 * @signature TVertexDescriptor getSource(stump);
 *
 * @param[in] stump Pointer to the EdgeStump to get the source of.
 *
 * @return TVertexDescriptor The vetex descriptor stored in stump.
 *
 * @section Remarks
 *
 * A source vertex is not required in an edge stump.  However, EdgeStump objects can be configured to contain a source
 * vertex, as in undirected graphs.
 */

template<typename TCargo, bool TList, bool TId, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TList, true, TId, TSpec> const>::Type
getSource(EdgeStump<TCargo, TList, true, TId, TSpec> const* es)
{
    SEQAN_CHECKPOINT
    return es->data_source;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TId, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TList, true, TId, TSpec> >::Type
getSource(EdgeStump<TCargo, TList, true, TId, TSpec>* es)
{
    SEQAN_CHECKPOINT
    return es->data_source;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TId, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TList, false, TId, TSpec> const>::Type
getSource(EdgeStump<TCargo, TList, false, TId, TSpec> const*)
{
    SEQAN_CHECKPOINT
    // Nop
    return 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TId, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TList, false, TId, TSpec> >::Type
getSource(EdgeStump<TCargo, TList, false, TId, TSpec>*)
{
    SEQAN_CHECKPOINT
    // Nop
    return 0;
}




//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn EdgeStump#assignNextT
 * @brief Assigns another EdgeStump to the next target pointer.
 *
 * @signature void assignNextT(es, es2);
 *
 * @param[in,out] es  Pointer to the EdgeStump.
 * @param[in]     es2 Pointer to the following EdgeStump.
 */

template<typename TCargo, bool TSource, bool TId, typename TSpec>
inline void
assignNextT(EdgeStump<TCargo, true, TSource, TId, TSpec>* es,
            EdgeStump<TCargo, true, TSource, TId, TSpec>* es2)
{
    SEQAN_CHECKPOINT
    es->data_nextT = es2;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn EdgeStump#nextT
 * @brief Accesses the next target pointer.
 *
 * @signature TEdgeStump nextT(es);
 *
 * @param[in] es Pointer to the EdgeStump.
 *
 * @return TEdgeStump Reference to the next target pointer.
 */

template<typename TCargo, bool TSource, bool TId, typename TSpec>
inline EdgeStump<TCargo, true, TSource, TId, TSpec>* &
nextT(EdgeStump<TCargo, true, TSource, TId, TSpec>* es)
{
    SEQAN_CHECKPOINT
    return es->data_nextT;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TSource, bool TId, typename TSpec>
inline EdgeStump<TCargo, true, TSource, TId, TSpec>* &
nextT(EdgeStump<TCargo, true, TSource, TId, TSpec> const* es)
{
    return es->data_nextT;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn EdgeStump#getNextT
 * @brief Get method for the next target pointer.
 *
 * @signature TEdgeStump getNextT(es);
 *
 * @param[in] es Pointer to the EdgeStump.
 *
 * @return TEdgeStump Reference to the next target pointer.
 */

template<typename TCargo, bool TSource, bool TId, typename TSpec>
inline EdgeStump<TCargo, true, TSource, TId, TSpec>*
getNextT(EdgeStump<TCargo, true, TSource, TId, TSpec>* es)
{
    SEQAN_CHECKPOINT
    return es->data_nextT;
}


//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TSource, bool TId, typename TSpec>
inline EdgeStump<TCargo, true, TSource, TId, TSpec>*
getNextT(EdgeStump<TCargo, true, TSource, TId, TSpec> const* es)
{
    return es->data_nextT;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn EdgeStump#assignNextS
 * @brief Assigns another EdgeStump to the next source pointer.
 *
 * @signature void assignNextS(es, es2);
 *
 * @param[in,out] es  Pointer to the EdgeStump.
 * @param[in]     es2 Pointer to the following EdgeStump.
 *
 * Edge Stumps can be configured to have no source.  In this case, there is no next source pointer.
 */

template<typename TCargo, bool TId, typename TSpec>
inline void
assignNextS(EdgeStump<TCargo, true, true, TId, TSpec>* es,
            EdgeStump<TCargo, true, true, TId, TSpec>* es2)
{
    SEQAN_CHECKPOINT
    es->data_nextS = es2;
}


//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TId, typename TSpec>
inline void
assignNextS(EdgeStump<TCargo, true, false, TId, TSpec>*,
            EdgeStump<TCargo, true, false, TId, TSpec>*)
{
    SEQAN_CHECKPOINT
    // Nop
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn EdgeStump#nextS
 * @brief Accesses the next source pointer.
 *
 * @signature TEdgeStump nextS(es);
 *
 * @param[in,out] es  Pointer to the EdgeStump.
 *
 * @return Reference to the next source pointer.
 *
 * Edge Stumps can be configured to have no source.  In this case, there is no next source pointer.
 */

template<typename TCargo, bool TId, typename TSpec>
inline EdgeStump<TCargo, true, true, TId, TSpec>* &
nextS(EdgeStump<TCargo, true, true, TId, TSpec>* es)
{
    SEQAN_CHECKPOINT
    return es->data_nextS;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TId, typename TSpec>
inline EdgeStump<TCargo, true, true, TId, TSpec>* &
nextS(EdgeStump<TCargo, true, true, TId, TSpec> const* es)
{
    return es->data_nextS;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TId, typename TSpec>
inline EdgeStump<TCargo, true, false, TId, TSpec>*
nextS(EdgeStump<TCargo, true, false, TId, TSpec>*)
{
    SEQAN_CHECKPOINT
    // Nop
    return 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TId, typename TSpec>
inline EdgeStump<TCargo, true, false, TId, TSpec>*
nextS(EdgeStump<TCargo, true, false, TId, TSpec> const*)
{
    // Nop
    return 0;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn EdgeStump#getNextS
 * @brief Accesses the next source pointer.
 *
 * @signature TEdgeStump getNextS(es);
 *
 * @param[in,out] es  Pointer to the EdgeStump.
 *
 * @return Reference to the next source pointer.
 *
 * Edge Stumps can be configured to have no source.  In this case, there is no next source pointer.
 */

template<typename TCargo, bool TId, typename TSpec>
inline EdgeStump<TCargo, true, true, TId, TSpec>*
getNextS(EdgeStump<TCargo, true, true, TId, TSpec> const* es)
{
    SEQAN_CHECKPOINT
    return es->data_nextS;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TId, typename TSpec>
inline EdgeStump<TCargo, true, false, TId, TSpec>*
getNextS(EdgeStump<TCargo, true, false, TId, TSpec> const*)
{
    SEQAN_CHECKPOINT
    // No source pointer
    return 0;
}




//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, typename TSpec, typename TId2>
void
_assignId(EdgeStump<TCargo, TList, TSource, true, TSpec>* es,
          TId2 const id)
{
    SEQAN_CHECKPOINT
    es->data_id = id;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, typename TSpec, typename TId2>
void
_assignId(EdgeStump<TCargo, TList, TSource, false, TSpec>*,
          TId2 const)
{
    // No id -> does nothing
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, typename TId2>
void
_assignId(EdgeStump<TCargo, TList, TSource, false, TreeTag>*,
          TId2 const)
{
    // For a tree do nothing, child id = tree id
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, typename TSpec>
inline typename Id<EdgeStump<TCargo, TList, TSource, true, TSpec> const>::Type
_getId(EdgeStump<TCargo, TList, TSource, true, TSpec> const* es)
{
    SEQAN_CHECKPOINT
    return es->data_id;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, typename TSpec>
inline typename Id<EdgeStump<TCargo, TList, TSource, true, TSpec> >::Type
_getId(EdgeStump<TCargo, TList, TSource, true, TSpec>* es)
{
    SEQAN_CHECKPOINT
    return es->data_id;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource>
inline typename Id<EdgeStump<TCargo, TList, TSource, false, TreeTag> const>::Type
_getId(EdgeStump<TCargo, TList, TSource, false, TreeTag> const* es)
{
    SEQAN_CHECKPOINT
    // Child id = edge id in a tree
    return es->data_target;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource>
inline typename Id<EdgeStump<TCargo, TList, TSource, false, TreeTag> >::Type
_getId(EdgeStump<TCargo, TList, TSource, false, TreeTag>* es)
{
    SEQAN_CHECKPOINT
    // Child id = edge id in a tree
    return es->data_target;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, typename TSpec>
inline typename Id<EdgeStump<TCargo, TList, TSource, false, TSpec> >::Type
_getId(EdgeStump<TCargo, TList, TSource, false, TSpec> const*)
{
    SEQAN_CHECKPOINT
    // No real id
    return 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, typename TSpec>
inline typename Id<EdgeStump<TCargo, TList, TSource, false, TSpec> >::Type
_getId(EdgeStump<TCargo, TList, TSource, false, TSpec>*)
{
    SEQAN_CHECKPOINT
    // No real id
    return 0;
}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
