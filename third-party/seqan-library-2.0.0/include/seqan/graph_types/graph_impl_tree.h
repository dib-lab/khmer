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

#ifndef SEQAN_HEADER_GRAPH_IMPL_TREE_H
#define SEQAN_HEADER_GRAPH_IMPL_TREE_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Graph - Tree
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/*!
 * @class Tree
 * @headerfile <seqan/graph_types.h>
 * @extends Graph
 * @brief A tree.
 *
 * A Tree has a distinct root and directed edges.  The source vertex of each edge is the parent vertex, the target
 * vertex of each edge is the child.  Trees provide fast access to child vertices and the parent.
 *
 * @signature template <[typename TCargo[, typename TSpec]]>
 *            class Graph<Tree<TCargo, TSpec> >;
 *
 * @tparam TCargo The cargo type that can be attached to the edges.  Default: <tt>void</tt>.
 * @tparam TSpec  Specializing type.  Default: <tt>Default</tt>.
 */

template<typename TCargo, typename TSpec>
class Graph<Tree<TCargo, TSpec> >
{
    public:
        typedef typename VertexIdHandler<Graph>::Type TVertexIdManager_;
        typedef typename VertexDescriptor<Graph>::Type TVertexDescriptor_;
        typedef typename EdgeType<Graph>::Type TEdgeStump_;
        typedef Allocator<SinglePool<sizeof(TEdgeStump_)> > TAllocator_;

        TVertexDescriptor_ data_root;
        String<TEdgeStump_*> data_vertex;            // Pointers to EdgeStumpT lists
        String<TVertexDescriptor_> data_parent;        // Map to the parents of each node
        TVertexIdManager_ data_id_managerV;
        TAllocator_ data_allocator;


//____________________________________________________________________________


        Graph() : data_root(getNil<TVertexDescriptor_>()) {
            SEQAN_CHECKPOINT
        }


        ~Graph() {
            SEQAN_CHECKPOINT
            clear(*this);
        }

        Graph(Graph const & _other) :
            data_root(getNil<TVertexDescriptor_>()),
            data_allocator(_other.data_allocator)
        {
            SEQAN_CHECKPOINT
            _copyGraph(_other, *this);
        }

        Graph const& operator = (Graph const & _other) {
            SEQAN_CHECKPOINT
            if (this == &_other) return *this;
            clear(*this);
            data_allocator = _other.data_allocator;
            _copyGraph(_other, *this);
            return *this;
        }
};


//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline String<typename EdgeType<Graph<Tree<TCargo, TSpec> > >::Type*>&
_getVertexString(Graph<Tree<TCargo, TSpec> > const& g) {
    SEQAN_CHECKPOINT
    typedef Graph<Tree<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    return const_cast<String<TEdgeStump*>&>(g.data_vertex);
}

/////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline IdManager<typename Id<Graph<Tree<TCargo, TSpec> > >::Type, Default> const &
_getVertexIdManager(Graph<Tree<TCargo, TSpec> > const& g) {
    return g.data_id_managerV;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline IdManager<typename Id<Graph<Tree<TCargo, TSpec> > >::Type, Default> &
_getVertexIdManager(Graph<Tree<TCargo, TSpec> >& g) {
    return g.data_id_managerV;
}

/////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline IdManager<typename Id<Graph<Tree<TCargo, TSpec> > >::Type, Default> const &
_getEdgeIdManager(Graph<Tree<TCargo, TSpec> > const& g) {
    return g.data_id_managerV;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline IdManager<typename Id<Graph<Tree<TCargo, TSpec> > >::Type, Default> &
_getEdgeIdManager(Graph<Tree<TCargo, TSpec> >& g) {
    return g.data_id_managerV;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline void
_rebuildParentMap(Graph<Tree<TCargo, TSpec> >& g)
{
    SEQAN_CHECKPOINT

    typedef Graph<Tree<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Iterator<String<TEdgeStump*>, Standard>::Type TIterConst;
    TIterConst it = begin(g.data_vertex, Standard());
    TIterConst itEnd = end(g.data_vertex, Standard());
    TVertexDescriptor pos = 0;
    for(;it!=itEnd; ++it, ++pos) {
        TVertexDescriptor parent = pos;
        TEdgeStump* current = *it;
        while(current!=0) {
            g.data_parent[getTarget(current)] = parent;
            current = getNextT(current);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline void
_copyGraph(Graph<Tree<TCargo, TSpec> > const& source,
           Graph<Tree<TCargo, TSpec> >& dest,
           bool transpose)
{
    SEQAN_CHECKPOINT
    typedef Graph<Tree<TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename Iterator<String<TEdgeStump*> const, Standard>::Type TIterConst;
    typedef typename Iterator<String<TEdgeStump*>, Standard>::Type TIter;
    clear(dest);
    resize(dest.data_vertex, length(source.data_vertex));
    resize(dest.data_parent, length(source.data_parent));
    TIter itInit = begin(dest.data_vertex, Standard());
    TIter itInitEnd = end(dest.data_vertex, Standard());
    for(;itInit != itInitEnd; ++itInit) *itInit = (TEdgeStump*) 0;
    TIterConst it = begin(source.data_vertex, Standard());
    TIterConst itEnd = end(source.data_vertex, Standard());
    TVertexDescriptor pos = 0;
    for(;it!=itEnd;++it, ++pos) {
        TEdgeStump* current = *it;
        TVertexDescriptor parentVertex = pos;
        while(current != (TEdgeStump*) 0) {
            TVertexDescriptor childVertex = getTarget(current);
            // Create missing vertices
            if (parentVertex > childVertex) _createVertices(dest,parentVertex);
            else _createVertices(dest,childVertex);
            // Add edge
            TEdgeDescriptor e;
            if (!transpose) e = addEdge(dest, parentVertex, childVertex);
            else e = addEdge(dest, childVertex, parentVertex);
            assignCargo(e, getCargo(current));
            current = getNextT(current);
        }
    }
    dest.data_id_managerV = source.data_id_managerV;
    dest.data_root = source.data_root;
}

template<typename TCargo, typename TSpec>
inline void
_copyGraph(Graph<Tree<TCargo, TSpec> > const& source,
           Graph<Tree<TCargo, TSpec> >& dest)
{
    _copyGraph(source, dest, false);
}

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline void
transpose(Graph<Tree<TCargo, TSpec> > const& source,
          Graph<Tree<TCargo, TSpec> >& dest)
{
    SEQAN_CHECKPOINT
    _copyGraph(source, dest, true);

}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline void
transpose(Graph<Tree<TCargo, TSpec> >& g)
{
    SEQAN_CHECKPOINT
    Graph<Tree<TCargo, TSpec> > dest;
    _copyGraph(g, dest, true);
    g = dest;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename Size<Graph<Tree<TCargo, TSpec> > >::Type
numEdges(Graph<Tree<TCargo, TSpec> > const& g)
{
    SEQAN_CHECKPOINT
    typedef Graph<Tree<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename Iterator<String<TEdgeStump*> const, Standard>::Type TIterConst;
    TSize count=0;
    TIterConst it = begin(g.data_vertex, Standard());
    TIterConst itEnd = end(g.data_vertex, Standard());
    for(;it!=itEnd;++it) {
        TEdgeStump* current = *it;
        while(current!=0) {
            ++count;
            current = getNextT(current);
        }
    }
    return count;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename Size<Graph<Tree<TCargo, TSpec> > >::Type
numVertices(Graph<Tree<TCargo, TSpec> > const& g)
{
    SEQAN_CHECKPOINT
    return idCount(g.data_id_managerV);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline bool
empty(Graph<Tree<TCargo, TSpec> > const& g)
{
    SEQAN_CHECKPOINT
    return (!idCount(g.data_id_managerV));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline void
clearEdges(Graph<Tree<TCargo, TSpec> >& g)
{
    SEQAN_CHECKPOINT

    typedef Graph<Tree<TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename Iterator<String<TEdgeStump*>, Standard>::Type TIter;
    TIter it = begin(g.data_vertex, Standard());
    TIter itEnd = end(g.data_vertex, Standard());
    TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
    TVertexDescriptor pos = 0;
    for(;it!=itEnd;++it, ++pos) {
        TEdgeStump* current = *it;
        if(current != (TEdgeStump*) 0) {
            g.data_parent[pos] = nilVertex;
            removeOutEdges(g, pos);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline void
clearVertices(Graph<Tree<TCargo, TSpec> >& g)
{
    SEQAN_CHECKPOINT
    clearEdges(g);
    releaseAll(g.data_id_managerV);
    clear(g.data_vertex);
    clear(g.data_parent);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline void
clear(Graph<Tree<TCargo, TSpec> >& g)
{
    SEQAN_CHECKPOINT
    clearVertices(g);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Tree<TCargo, TSpec> > >::Type
outDegree(Graph<Tree<TCargo, TSpec> > const& g,
          TVertexDescriptor const vertex)
{
    SEQAN_CHECKPOINT
    SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex));

    typedef Graph<Tree<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename Size<TGraph>::Type TSize;
    TSize count=0;
    TEdgeStump* current = g.data_vertex[vertex];
    while(current!=0) {
        current = getNextT(current);
        ++count;
    }
    return count;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Tree<TCargo, TSpec> > >::Type
inDegree(Graph<Tree<TCargo, TSpec> > const& g,
         TVertexDescriptor const vertex)
{
    SEQAN_CHECKPOINT
    SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex));

    typedef Graph<Tree<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename Iterator<String<TEdgeStump*> const, Standard>::Type TIterConst;
    TSize count=0;
    TIterConst it = begin(g.data_vertex, Standard());
    TIterConst itEnd = end(g.data_vertex, Standard());
    for(;it!=itEnd;++it) {
        TEdgeStump* current = *it;
        while(current!=0) {
            if ( (TVertexDescriptor) getTarget(current)==vertex) ++count;
            current = getNextT(current);
        }
    }
    return count;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Tree<TCargo, TSpec> > >::Type
degree(Graph<Tree<TCargo, TSpec> > const& g,
       TVertexDescriptor const vertex)
{
    SEQAN_CHECKPOINT
    return (inDegree(g,vertex)+outDegree(g,vertex));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<Graph<Tree<TCargo, TSpec> > >::Type
addVertex(Graph<Tree<TCargo, TSpec> >& g)
{
    SEQAN_CHECKPOINT
    typedef Graph<Tree<TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
    TVertexDescriptor vd;
    if (empty(g)) g.data_root = vd = obtainId(g.data_id_managerV);
    else vd = obtainId(g.data_id_managerV);
    if (vd == length(g.data_vertex)) {
        appendValue(g.data_vertex, (TEdgeStump*) 0);
        resize(g.data_parent, vd + 1, nilVertex, Generous());
    } else {
        g.data_vertex[vd] = (TEdgeStump*) 0;
        g.data_parent[vd] = nilVertex;
    }
    return vd;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
removeVertex(Graph<Tree<TCargo, TSpec> >& g,
             TVertexDescriptor const v)
{
    SEQAN_CHECKPOINT
    SEQAN_ASSERT(idInUse(g.data_id_managerV, v));

    typedef Graph<Tree<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
    TEdgeStump* current = getValue(g.data_vertex, v);
    while(current!=0) {
        g.data_parent[childVertex(g, current)] = nilVertex;
        current = getNextT(current);
    }
    g.data_parent[v] = nilVertex;
    removeOutEdges(g,v); // Remove all outgoing edges
    removeInEdges(g,v); // Remove all incoming edges
    releaseId(g.data_id_managerV, v); // Release id
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename EdgeDescriptor<Graph<Tree<TCargo, TSpec> > >::Type
addEdge(Graph<Tree<TCargo, TSpec> >& g,
        TVertexDescriptor const parent,
        TVertexDescriptor const child)
{
    SEQAN_ASSERT(idInUse(g.data_id_managerV, parent));
    SEQAN_ASSERT(idInUse(g.data_id_managerV, child));
    SEQAN_ASSERT(findEdge(g, parent, child) == 0); // No multi-graphs as trees!!!

    typedef Graph<Tree<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;

    TEdgeStump* edge_ptr;
    allocate(g.data_allocator, edge_ptr, 1);
    valueConstruct(edge_ptr);
    assignTarget(edge_ptr, child);
    g.data_parent[child] = parent;
    assignNextT(edge_ptr, (TEdgeStump*) 0);
    if (g.data_vertex[parent]!=0) assignNextT(edge_ptr, g.data_vertex[parent]);
    g.data_vertex[parent]=edge_ptr;
    return edge_ptr;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename EdgeDescriptor<Graph<Tree<TCargo, TSpec> > >::Type
addEdge(Graph<Tree<TCargo, TSpec> >& g,
        TVertexDescriptor const parent,
        TVertexDescriptor const child,
        TCargo const cargo)
{
    SEQAN_CHECKPOINT
    typedef Graph<Tree<TCargo, TSpec> > TGraph;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    TEdgeDescriptor e = addEdge(g,parent,child);
    assignCargo(e,cargo);
    return e;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
removeEdge(Graph<Tree<TCargo, TSpec> >& g,
           TVertexDescriptor const parent,
           TVertexDescriptor const child)
{
    SEQAN_CHECKPOINT
    SEQAN_ASSERT(idInUse(g.data_id_managerV, parent));
    SEQAN_ASSERT(idInUse(g.data_id_managerV, child));

    typedef Graph<Tree<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;

    // Find edge and predecessor
    TEdgeStump* pred = 0;
    TEdgeStump* current = g.data_vertex[parent];
    while(current != (TEdgeStump*) 0) {
        if ( (TVertexDescriptor) getTarget(current) == child) break;
        pred = current;
        current = getNextT(current);
    }

    // Not found?
    if (current == (TEdgeStump*) 0) return;
    g.data_parent[child] = getNil<TVertexDescriptor>();

    // Relink the next pointer of predecessor
    if (pred != (TEdgeStump*) 0) assignNextT(pred, getNextT(current));
    else g.data_vertex[parent] = getNextT(current);

    // Deallocate
    valueDestruct(current);
    deallocate(g.data_allocator, current, 1);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline void
removeEdge(Graph<Tree<TCargo, TSpec> >& g,
           TEdgeDescriptor const edge)
{
    SEQAN_CHECKPOINT
    removeEdge(g, parentVertex(g,edge), childVertex(g,edge));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
removeOutEdges(Graph<Tree<TCargo, TSpec> >& g,
               TVertexDescriptor const v)
{
    SEQAN_CHECKPOINT
    SEQAN_ASSERT(idInUse(g.data_id_managerV, v));

    typedef Graph<Tree<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    while(g.data_vertex[v] != (TEdgeStump*) 0) {
        TVertexDescriptor target = targetVertex(g,g.data_vertex[v]);
        removeEdge(g,v,target);
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
removeInEdges(Graph<Tree<TCargo, TSpec> >& g,
              TVertexDescriptor const v)
{
    SEQAN_CHECKPOINT
    SEQAN_ASSERT(idInUse(g.data_id_managerV, v));

    typedef Graph<Tree<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename Iterator<String<TEdgeStump*>, Standard>::Type TIter;
    TIter it = begin(g.data_vertex, Standard());
    TIter itEnd = end(g.data_vertex, Standard());
    TVertexDescriptor pos = 0;
    for(;it!=itEnd;++it, ++pos) {
        if (!idInUse(g.data_id_managerV, pos)) continue;
        TEdgeStump* current = *it;
        TVertexDescriptor const sourceVertex = pos;
        while(current!=0) {
            if ( (TVertexDescriptor) current->data_target==v) {
                removeEdge(g, sourceVertex, v);
                current = g.data_vertex[sourceVertex];
            } else {
                current = getNextT(current);
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Tree<TCargo, TSpec> > >::Type
targetVertex(Graph<Tree<TCargo, TSpec> > const&,
             TEdgeDescriptor const edge)
{
    SEQAN_CHECKPOINT
    return getTarget(edge);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Tree<TCargo, TSpec> > >::Type
sourceVertex(Graph<Tree<TCargo, TSpec> > const& g,
             TEdgeDescriptor const edge)
{
    SEQAN_CHECKPOINT
    typedef Graph<Tree<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Iterator<String<TEdgeStump*> const, Standard>::Type TIterConst;
    TIterConst it = begin(g.data_vertex, Standard());
    TIterConst itEnd = end(g.data_vertex, Standard());
    TVertexDescriptor pos = 0;
    for(;it!=itEnd;++it, ++pos) {
        TEdgeDescriptor current = *it;
        while(current!=(TEdgeDescriptor) 0) {
            if (current == edge) return pos;
            current=getNextT(current);
        }
    }
    return 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TMatrix>
inline void
getAdjacencyMatrix(Graph<Tree<TCargo, TSpec> > const& g,
                   TMatrix& mat)
{
    SEQAN_CHECKPOINT
    typedef Graph<Tree<TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename Size<TMatrix>::Type TSize;
    TSize len = getIdUpperBound(g.data_id_managerV);
    resize(mat, len*len, 0);

    typedef typename Iterator<String<TEdgeStump*> const, Standard>::Type TIterConst;
    TIterConst it = begin(g.data_vertex, Standard());
    TIterConst itEnd = end(g.data_vertex, Standard());
    TVertexDescriptor pos = 0;
    for(;it!=itEnd;++it,++pos) {
        TVertexDescriptor parentV = pos;
        TEdgeStump* current = *it;
        while(current!=0) {
            TVertexDescriptor childV = getTarget(current);
            ++mat[parentV*len+childV];
            current=getNextT(current);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename EdgeDescriptor<Graph<Tree<TCargo, TSpec> > >::Type
findEdge(Graph<Tree<TCargo, TSpec> >& g,
         TVertexDescriptor const v,
         TVertexDescriptor const w)
{
    SEQAN_CHECKPOINT
    SEQAN_ASSERT(idInUse(g.data_id_managerV, v));
    SEQAN_ASSERT(idInUse(g.data_id_managerV, w));

    typedef Graph<Tree<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;

    if ((TVertexDescriptor)g.data_parent[w] == v) {
        TEdgeStump* current = g.data_vertex[v];
        while((TEdgeStump*) current != 0) {
            if ((TVertexDescriptor) getTarget(current) == w) return current;
            current = getNextT(current);
        }
    } else if ((TVertexDescriptor) g.data_parent[v] == w) {
        TEdgeStump* current = g.data_vertex[w];
        while((TEdgeStump*) current != 0) {
            if ((TVertexDescriptor) getTarget(current) == v) return current;
            current = getNextT(current);
        }
    }
    return 0;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
write(TFile & target,
      Graph<Tree<TCargo, TSpec> > const& g)
{
//IOREV _nodoc_
    SEQAN_CHECKPOINT
    typedef Graph<Tree<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Iterator<String<TEdgeStump*> const, Standard>::Type TIterConst;
    TIterConst it = begin(g.data_vertex, Standard());
    TIterConst itEnd = end(g.data_vertex, Standard());
    TVertexDescriptor pos = 0;
    write(target, "Adjacency list:\n");
    for(;it!=itEnd;++it, ++pos) {
        if (!idInUse(_getVertexIdManager(g),pos)) continue;
        TEdgeStump* current = *it;
        appendNumber(target, (int)pos);
        write(target, " -> ");
        while(current!=0) {
            appendNumber(target, (int)getTarget(current));
            writeValue(target, ',');
            current=getNextT(current);
        }
        writeValue(target, '\n');
    }
    it = begin(g.data_vertex, Standard());
    pos = 0;
    write(target, "Edge list:\n");
    for(;it!=itEnd;++it, ++pos) {
        TEdgeStump* current = getValue(it);
        while(current!=0) {
            write(target, "Source: ");
            appendNumber(target, (int)pos);
            writeValue(target, ',');
            write(target, "Target: ");
            appendNumber(target, (int)getTarget(current));
            writeValue(target, ' ');
            write(target, "(Id: ");
            appendNumber(target, (int)_getId(current));
            writeValue(target, ')');
            writeValue(target, '\n');
            current=getNextT(current);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Tree#assignRoot
 * @brief Assign a new root vertex to the tree.
 *
 * @signature void assignRoot(t, v);
 *
 * @param[in,out] t The Tree to assign the root for.
 * @param[in]     v A vertex descriptor.
 */

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
assignRoot(Graph<Tree<TCargo, TSpec> >& g,
           TVertexDescriptor const vertex)
{
    SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex));

    g.data_root = vertex;
    g.data_parent[vertex] = getNil<TVertexDescriptor>();
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Tree#createRoot
 * @brief Creates the root for the tree.
 *
 * @signature void createRoot(t);
 *
 * @param[in,out] t The Tree to create the root for.
 */

template<typename TCargo, typename TSpec>
inline void
createRoot(Graph<Tree<TCargo, TSpec> >& g)
{
    SEQAN_CHECKPOINT
    g.data_root = addVertex(g);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Tree#root
 * @brief Gets reference to the root of the tree.
 *
 * @signature TVertexDescriptor root(t);
 *
 * @param[in] t The Tree to query for its root.
 *
 * @return TVertexDescriptor Reference to the root's vertex descriptor.
 */

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<Graph<Tree<TCargo, TSpec> > >::Type&
root(Graph<Tree<TCargo, TSpec> >& g)
{
    SEQAN_CHECKPOINT
    return g.data_root;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Tree#getRoot
 * @brief Gets the root of the tree.
 *
 * @signature TVertexDescriptor getRoot(t);
 *
 * @param[in] t The Tree to query for its root.
 *
 * @return TVertexDescriptor The root's vertex descriptor.
 */

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<Graph<Tree<TCargo, TSpec> > >::Type
getRoot(Graph<Tree<TCargo, TSpec> > const& g)
{
    SEQAN_CHECKPOINT
    return g.data_root;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Tree#isRoot
 * @brief Tests whether a given vertex is the root or not.
 *
 * @signature bool isRoot(t, v);
 *
 * @param[in] a The Tree to query.
 * @param[in] t The descriptor of the vertex to query.
 *
 * @return bool true if v is the root and false otherwise.
 */

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline bool
isRoot(Graph<Tree<TCargo, TSpec> > const& g,
       TVertexDescriptor v)
{
    SEQAN_CHECKPOINT
    return ( (TVertexDescriptor) g.data_root == v);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Tree#isLeaf
 * @brief Tests whether a given vertex is a leaf or not.
 *
 * @signature bool isLeaf(g, v);
 *
 * @param[in] g The Tree to query.
 * @param[in] v The descriptor of the vertex to query for.
 */

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline bool
isLeaf(Graph<Tree<TCargo, TSpec> > const& g,
       TVertexDescriptor v)
{
    SEQAN_CHECKPOINT
    SEQAN_ASSERT(idInUse(g.data_id_managerV, v));

    typedef Graph<Tree<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    return (g.data_vertex[v] ==  (TEdgeStump*) 0);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Tree#numTreeEdges
 * @brief Number of tree edges.
 *
 * @signature TSize numTreeEdges(g);
 *
 * @param[in] g The tree to query.
 *
 * @return TSize The number of tree edges.  Faster than @link Graph#numEdges @endlink for trees.
 */

template<typename TCargo, typename TSpec>
inline typename Size<Graph<Tree<TCargo, TSpec> > >::Type
numTreeEdges(Graph<Tree<TCargo, TSpec> > const& g)
{
    SEQAN_CHECKPOINT
    if (empty(g)) return 0;
    else return numVertices(g) - 1;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Tree#numChildren
 * @brief Number of children of a given tree vertex.
 *
 * @signature TSize numChildren(g, v);
 *
 * @param[in] g The Tree to query.
 * @param[in] v The descriptor of the vertex to count the children of.
 *
 * @return TSize The number of children.
 */

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Tree<TCargo, TSpec> > >::Type
numChildren(Graph<Tree<TCargo, TSpec> > const& g,
            TVertexDescriptor const vertex)
{
    SEQAN_CHECKPOINT

    return outDegree(g, vertex);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Tree#addChild
 * @brief Adds a new child vertex to a parent vertex, optionally with an edge cargo.
 *
 * @signature TVertexDescriptor(g, parent[, cargo]);
 *
 * @param[in,out] g      The Tree to append the vertex and edge to.
 * @param[in]     parent Descriptor of the parent vertex.
 * @param[in]     cargo  The cargo to attach to the edge.
 *
 * @return TVertexDescriptor Vertex descriptor of the added vertex.
 */

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename VertexDescriptor<Graph<Tree<TCargo, TSpec> > >::Type
addChild(Graph<Tree<TCargo, TSpec> >& g,
         TVertexDescriptor parent)
{
    SEQAN_CHECKPOINT
    SEQAN_ASSERT(idInUse(g.data_id_managerV, parent));
    TVertexDescriptor child = addVertex(g);
    addEdge(g,parent,child);
    return child;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename VertexDescriptor<Graph<Tree<TCargo, TSpec> > >::Type
addChild(Graph<Tree<TCargo, TSpec> >& g,
         TVertexDescriptor const parent,
         TCargo const cargo)
{
    SEQAN_CHECKPOINT
    SEQAN_ASSERT(idInUse(g.data_id_managerV, parent));
    TVertexDescriptor child = addVertex(g);
    addEdge(g,parent,child,cargo);
    return child;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Tree#removeChild
 * @brief Removes a child from the tree given a parent.
 *
 * @signature void removeChidl(g, parent, child);
 *
 * @param[in,out] g      The Tree to remove the vertex from.
 * @param[in]     parent The descriptor of the parent vertex.
 * @param[in]     child  The descriptor of the child vertex.
 */

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
removeChild(Graph<Tree<TCargo, TSpec> >& g,
            TVertexDescriptor const parent,
            TVertexDescriptor const child)
{
    SEQAN_CHECKPOINT
    if (!isLeaf(g,child)) removeAllChildren(g,child);
    removeEdge(g,parent, child);  // Parent map is cleared in removeEdge
    releaseId(g.data_id_managerV, child); // Release id
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Tree#removeAllChildren
 * @brief Removes all children from the tree given a parent.
 *
 * @signature void removeAllChildren(g, parent);
 *
 * @param[in,out] g      The Tree to modify.
 * @param[in]     parent Descriptor of the vertex to remove children of.
 */

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
removeAllChildren(Graph<Tree<TCargo, TSpec> >& g,
                  TVertexDescriptor const parent)
{
    SEQAN_CHECKPOINT
    SEQAN_ASSERT(idInUse(g.data_id_managerV, parent));
    typedef Graph<Tree<TCargo, TSpec> > TGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;

    while(g.data_vertex[parent] != (TEdgeStump*) 0) {
        TVertexDescriptor child = childVertex(g,g.data_vertex[parent]);
        if (!isLeaf(g,child)) removeAllChildren(g,child);
        removeEdge(g,parent, child);  // Parent map is cleared in removeEdge
        releaseId(g.data_id_managerV, child); // Release id
    }
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Tree#childVertex
 * @brief Returns the child vertex of an edge.
 *
 * @signature TVertexDescriptor childVertex(g, e);
 *
 * @param[in] g The Tree to query.
 * @param[in] e The edge descriptor of the edge to query.
 *
 * @return TVertexDescriptor The descriptor of the child vertex.
 */

template<typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Tree<TCargo, TSpec> > >::Type
childVertex(Graph<Tree<TCargo, TSpec> > const&,
            TEdgeDescriptor const edge)
{
    SEQAN_CHECKPOINT
    return getTarget(edge);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Tree#parentVertex
 * @brief Returns the parent vertex of an edge.
 *
 * @signature TVertexDescriptor childVertex(g, e);
 * @signature TVertexDescriptor childVertex(g, v);
 *
 * @param[in] g The Tree to query.
 * @param[in] e The edge descriptor of the edge to query.
 * @param[in] v The descriptor of the vertex to query.
 *
 * @return TVertexDescriptor The descriptor of the child vertex.
 */

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<Graph<Tree<TCargo, TSpec> > >::Type
parentVertex(Graph<Tree<TCargo, TSpec> > const& g,
             typename EdgeDescriptor<Graph<Tree<TCargo, TSpec> > >::Type const edge)
{
    SEQAN_CHECKPOINT
    return g.data_parent[getTarget(edge)];
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<Graph<Tree<TCargo, TSpec> > >::Type
parentVertex(Graph<Tree<TCargo, TSpec> > const& g,
             typename VertexDescriptor<Graph<Tree<TCargo, TSpec> > >::Type const v)
{
    SEQAN_CHECKPOINT
    return g.data_parent[v];
}

///////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Tree#collectLeaves
 * @brief Returns all leaves underneath a given vertex.
 *
 * @signature void collectLeaves(g, v, group);
 *
 * @param[in]  g     The Tree to collect the leaves for.
 * @param[in]  v     The root of the subtree to query the leaves of.
 * @param[out] group A @link String @endlink of leaf vertex descriptors.
 */

template<typename TCargo, typename TSpec, typename TVertexDescriptor, typename TGroup>
inline void
collectLeaves(Graph<Tree<TCargo, TSpec> > const& g,
              TVertexDescriptor const root,
              TGroup& group)
{
    typedef Graph<Tree<TCargo, TSpec> > TGraph;
    typedef typename Iterator<TGraph, AdjacencyIterator>::Type TAdjacencyIterator;

    if (isLeaf(g, root)) appendValue(group, root, Generous());
    else {
        TAdjacencyIterator adjIt(g, root);
        for(;!atEnd(adjIt);goNext(adjIt)) {
            collectLeaves(g, *adjIt, group);
        }
    }
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
