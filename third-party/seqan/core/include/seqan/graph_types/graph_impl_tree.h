// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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

/**
.Spec.Tree:
..cat:Graph
..general:Class.Graph
..summary:A Tree has a distinct root and directed edges. The source vertex of each edge is the parent vertex, 
the target vertex of each edge is the child. Trees provide fast access to child vertices and the parent.
..description:
...image:treeGraph|A tree, where $0$ is the root vertex.
..signature:Graph<Tree<TCargo, TSpec> >
..param.TCargo:The cargo type that can be attached to the edges.
...metafunction:Metafunction.Cargo
...remarks:Use @Metafunction.Cargo@ to get the cargo type of the tree.
...default:$void$
..param.TSpec:The specializing type for the graph.
...metafunction:Metafunction.Spec
...default:$Default$, see @Tag.Default@.
..include:seqan/graph_types.h
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
		String<TEdgeStump_*> data_vertex;			// Pointers to EdgeStumpT lists
		String<TVertexDescriptor_> data_parent;		// Map to the parents of each node
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

template <typename TFile, typename TCargo, typename TSpec, typename TIDString>
inline void
write(TFile & target,
	  Graph<Tree<TCargo, TSpec> > const& g,
	  TIDString const &,
	  Raw)
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
	streamPut(target,"Adjacency list:\n");
	for(;it!=itEnd;++it, ++pos) {
		if (!idInUse(_getVertexIdManager(g),pos)) continue;
		TEdgeStump* current = *it;
		streamPut(target, (int)pos);
		streamPut(target," -> ");
		while(current!=0) {
			streamPut(target, (int)getTarget(current));
			streamPut(target, ',');
			current=getNextT(current);
		}
		streamPut(target, '\n');
	}
	it = begin(g.data_vertex, Standard());
	pos = 0;
	streamPut(target,"Edge list:\n");
	for(;it!=itEnd;++it, ++pos) {
		TEdgeStump* current = getValue(it);
		while(current!=0) {
			streamPut(target,"Source: ");
			streamPut(target, (int)pos);		
			streamPut(target, ',');
			streamPut(target,"Target: ");
			streamPut(target, (int)getTarget(current));
			streamPut(target, ' ');
			streamPut(target,"(Id: ");
			streamPut(target, (int)_getId(current));
			streamPut(target, ')');
			streamPut(target, '\n');
			current=getNextT(current);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

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

template<typename TCargo, typename TSpec>
inline void
createRoot(Graph<Tree<TCargo, TSpec> >& g)
{
	SEQAN_CHECKPOINT
	g.data_root = addVertex(g);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<Graph<Tree<TCargo, TSpec> > >::Type&
root(Graph<Tree<TCargo, TSpec> >& g)
{
	SEQAN_CHECKPOINT
	return g.data_root;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<Graph<Tree<TCargo, TSpec> > >::Type 
getRoot(Graph<Tree<TCargo, TSpec> > const& g)
{
	SEQAN_CHECKPOINT
	return g.data_root;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline bool
isRoot(Graph<Tree<TCargo, TSpec> > const& g,
	   TVertexDescriptor v)
{
	SEQAN_CHECKPOINT
	return ( (TVertexDescriptor) g.data_root == v);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.Graph#isLeaf
..class:Spec.Tree
..cat:Graph
..summary:Tests whether a given vertex is a leaf or not.
..signature:isLeaf(g, v)
..param.g:A tree.
...type:Spec.Tree
..param.v:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:True if vertex is a leaf.
..include:seqan/graph_types.h
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

/**
.Function.numTreeEdges
..class:Spec.Tree
..cat:Graph
..summary:Number of tree edges.
..signature:numTreeEdges(g)
..param.g:A tree.
...type:Spec.Tree
..returns:Number of tree edges. Faster than numEdges for trees.
..include:seqan/graph_types.h
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

/**
.Function.numChildren:
..class:Spec.Tree
..cat:Graph
..summary:Number of children of a given tree vertex.
..signature:numChildren(g, v)
..param.g:A tree.
...type:Spec.Tree
..param.v:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:Number of children
..include:seqan/graph_types.h
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

/**
.Function.addChild
..class:Spec.Tree
..cat:Graph
..summary:Adds a new child vertex to a parent vertex.
Optionally a cargo can be attached to the parent-child edge.
..signature:addChild(g, parent [, cargo])
..param.g:A tree.
...type:Spec.Tree
..param.parent:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..param.cargo:A cargo object.
...type:Metafunction.Cargo
..returns:A new vertex descriptor.
...type:Metafunction.VertexDescriptor
..see:Function.removeChild
..see:Function.removeAllChildren
..include:seqan/graph_types.h
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

/**
.Function.removeChild
..class:Spec.Tree
..cat:Graph
..summary:Removes a child from the tree given a parent.
..signature:removeChild(g, parent, child)
..param.g:A tree.
...type:Spec.Tree
..param.parent:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..param.child:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:void
..see:Function.addChild
..see:Function.removeAllChildren
..include:seqan/graph_types.h
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

/**
.Function.removeAllChildren
..class:Spec.Tree
..cat:Graph
..summary:Removes all children from the tree given a parent.
..signature:removeChild(g, parent)
..param.g:A tree.
...type:Spec.Tree
..param.parent:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:void
..see:Function.addChild
..see:Function.removeChild
..include:seqan/graph_types.h
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

/**
.Function.childVertex
..class:Spec.Tree
..cat:Graph
..summary:Returns the child vertex of an edge.
..signature:childVertex(g, e)
..param.g:A tree.
...type:Spec.Tree
..param.e:An edge descriptor.
...type:Metafunction.EdgeDescriptor
..returns:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..see:Function.parentVertex
..include:seqan/graph_types.h
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

/**
.Function.parentVertex
..class:Spec.Tree
..cat:Graph
..summary:Returns the parent vertex of an edge or vertex.
..signature:parentVertex(g, e)
..signature:parentVertex(g, v)
..param.g:A tree.
...type:Spec.Tree
..param.e:An edge descriptor.
...type:Metafunction.EdgeDescriptor
..param.v:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..see:Function.parentVertex
..include:seqan/graph_types.h
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


/**
.Function.collectLeaves
..class:Spec.Tree
..cat:Graph
..summary:Returns all leaves underneath a given vertex.
..signature:collectLeaves(g, subtree_root, group)
..param.g:A tree.
...type:Spec.Tree
..param.subtree_root:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..param.group:A group of vertex descriptors
...type:Class.String
..include:seqan/graph_types.h
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
