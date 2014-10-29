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

#ifndef SEQAN_HEADER_GRAPH_ITERATOR_OUTEDGE_H
#define SEQAN_HEADER_GRAPH_ITERATOR_OUTEDGE_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Graph OutEdgeIterator
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Out-Edge Iterator:
..cat:Graph
..summary:Out-edge iterator for @Class.Graph@.
..signature:Iterator<TGraph, OutEdgeIterator>
..param.TGraph:A graph.
...type:Class.Graph
..general:Class.Iter
..see:Spec.Vertex Iterator
..see:Spec.Edge Iterator
..see:Spec.Adjacency Iterator
..see:Spec.Bfs Iterator
..see:Spec.Dfs Preorder Iterator
..include:seqan/graph_types.h

.Memfunc.Out-Edge Iterator#Iter
..summary:Constructor
..class:Spec.Out-Edge Iterator
..signature:Iter()
..signature:Iter(graph, v)
..signature:Iter(iter)
..param.graph:The graph to iterate edges of.
...type:Class.Graph
..param.v:Vertex descriptor of vertex to iterate out-edges from.
...type:Metafunction.VertexDescriptor
..param.iter:Other iterator to copy-construct from.
...type:Spec.Out-Edge Iterator
*/


//////////////////////////////////////////////////////////////////////////////
// InternalOutEdgeIterator for Directed Graph
//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
class Iter<Graph<Directed<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > > 
{
public:
	typedef Graph<Directed<TCargo, TGraphSpec> > TGraph_;
	typedef typename EdgeDescriptor<TGraph_>::Type TEdgeDescriptor_;
	typedef typename VertexDescriptor<TGraph_>::Type TVertexDescriptor_;
	TGraph_ const* data_host;
	TVertexDescriptor_ data_source;
	TEdgeDescriptor_ data_edge;

	Iter()	
	{
		SEQAN_CHECKPOINT
	}
	
	Iter(TGraph_ const& _graph, TVertexDescriptor_ const v) : 
		data_host(&_graph),
		data_source(v),
		data_edge(getValue(_graph.data_vertex,v))
	{
		SEQAN_CHECKPOINT
	}
	
	Iter(Iter const& _iter) : 
		data_host(_iter.data_host),
		data_source(_iter.data_source),
		data_edge(_iter.data_edge)
	{
		SEQAN_CHECKPOINT
	}

	~Iter() {
		SEQAN_CHECKPOINT
	}

	Iter const&	operator = (Iter const & _other) {
		SEQAN_CHECKPOINT
		if (this == &_other) return *this;
		data_host = _other.data_host;
		data_source = _other.data_source;
		data_edge = _other.data_edge;
		return *this;
	}
//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
// Graph InternalOutEdgeIterator - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec>
struct Iterator<Graph<Directed<TCargo, TGraphSpec> >, OutEdgeIterator>
{	
	typedef Iter<Graph<Directed<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<OutEdgeIterator> > > Type;
};

template<typename TCargo, typename TGraphSpec>
struct Iterator<Graph<Directed<TCargo, TGraphSpec> > const, OutEdgeIterator> :
    public Iterator<Graph<Directed<TCargo, TGraphSpec> >, OutEdgeIterator> {};


//////////////////////////////////////////////////////////////////////////////
// InternalOutEdgeIterator for Tree
//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
class Iter<Graph<Tree<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > > 
{
public:
	typedef Graph<Tree<TCargo, TGraphSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TGraph const* data_host;
	TVertexDescriptor data_source;
	TEdgeDescriptor data_edge;

	Iter()	
	{
		SEQAN_CHECKPOINT
	}
	
    Iter(TGraph const& _graph, TVertexDescriptor const v) :
		data_host(&_graph),
		data_source(v),
		data_edge(getValue(_graph.data_vertex,v))
	{
		SEQAN_CHECKPOINT
	}
	
	Iter(Iter const& _iter) :
		data_host(_iter.data_host),
		data_source(_iter.data_source),
		data_edge(_iter.data_edge)
	{
		SEQAN_CHECKPOINT
	}

	~Iter() {
		SEQAN_CHECKPOINT
	}

	Iter const&	operator = (Iter const & _other) {
		SEQAN_CHECKPOINT
		if (this == &_other) return *this;
		data_host = _other.data_host;
		data_source = _other.data_source;
		data_edge = _other.data_edge;
		return *this;
	}
//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
// Graph InternalOutEdgeIterator - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec>
struct Iterator<Graph<Tree<TCargo, TGraphSpec> >, OutEdgeIterator>
{	
	typedef Iter<Graph<Tree<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<OutEdgeIterator> > > Type;
};

template<typename TCargo, typename TGraphSpec>
struct Iterator<Graph<Tree<TCargo, TGraphSpec> > const, OutEdgeIterator> :
    public Iterator<Graph<Tree<TCargo, TGraphSpec> >, OutEdgeIterator> {};


//////////////////////////////////////////////////////////////////////////////
// InternalOutEdgeIterator for Undirected graph
//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
class Iter<Graph<Undirected<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > > 
{
public:
	typedef Graph<Undirected<TCargo, TGraphSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TGraph const* data_host;
	TVertexDescriptor data_source;
	TEdgeDescriptor data_edge;

	Iter()	
	{
		SEQAN_CHECKPOINT
	}
	
	Iter(TGraph const& _graph, TVertexDescriptor const v) : 
		data_host(&_graph),
		data_source(v),
		data_edge(getValue(_graph.data_vertex,v))
	{
		SEQAN_CHECKPOINT
	}
	
	Iter(Iter const& _iter) : 
		data_host(_iter.data_host),
		data_source(_iter.data_source),
		data_edge(_iter.data_edge)
	{
		SEQAN_CHECKPOINT
	}

	~Iter() {
		SEQAN_CHECKPOINT
	}

	Iter const&	operator = (Iter const & _other) {
		SEQAN_CHECKPOINT
		if (this == &_other) return *this;
		data_host = _other.data_host;
		data_source = _other.data_source;
		data_edge = _other.data_edge;
		return *this;
	}
//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
// Graph InternalOutEdgeIterator - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec>
struct Iterator<Graph<Undirected<TCargo, TGraphSpec> >, OutEdgeIterator>
{	
	typedef Iter<Graph<Undirected<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<OutEdgeIterator> > > Type;
};

template<typename TCargo,typename TGraphSpec>
struct Iterator<Graph<Undirected<TCargo, TGraphSpec> > const, OutEdgeIterator> :
    public Iterator<Graph<Undirected<TCargo, TGraphSpec> >, OutEdgeIterator> {};



//////////////////////////////////////////////////////////////////////////////
// Graph InternalOutEdgeIterator for Automaton
//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TIteratorSpec>
class Iter<Graph<Automaton<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TIteratorSpec> > > 
{
public:
	typedef Graph<Automaton<TAlphabet, TCargo, TGraphSpec> > TGraph_;
	typedef typename EdgeDescriptor<TGraph_>::Type TEdgeDescriptor_;
	typedef typename VertexDescriptor<TGraph_>::Type TVertexDescriptor_;
	typedef typename Size<TAlphabet>::Type TSize_;
	TGraph_ const* data_host;
	TVertexDescriptor_ data_source;
	TSize_ data_pos;
	TSize_ data_begin;
	TSize_ data_end;

	Iter()	
	{
		SEQAN_CHECKPOINT
	}

	Iter(TGraph_ const& _graph, TVertexDescriptor_ const v) : 
		data_host(&_graph),
		data_source(v)
	{
		SEQAN_CHECKPOINT
		TVertexDescriptor_ nilVal = getNil<TVertexDescriptor_>();
		TSize_ table_length = ValueSize<TAlphabet>::VALUE;
		TSize_ pos = 0;
		while (	(pos < table_length) &&
				(_graph.data_vertex[v].data_edge[pos].data_target == nilVal))
		{
				++pos;
		}
		data_pos = pos;
		data_begin = pos;
		data_end = table_length;
	}

	Iter(Iter const& _iter) : 
		data_host(_iter.data_host),
		data_source(_iter.data_source),
		data_pos(_iter.data_pos),
		data_begin(_iter.data_begin),
		data_end(_iter.data_end)
	{
		SEQAN_CHECKPOINT
	}

	~Iter() 
	{
		SEQAN_CHECKPOINT
	}

	Iter const&	operator = (Iter const & _other) 
	{
		SEQAN_CHECKPOINT
		if (this == &_other) return *this;
		data_host = _other.data_host;
		data_source = _other.data_source;
		data_pos = _other.data_pos;
		data_begin = _other.data_begin;
		data_end = _other.data_end;
		return *this;
	}
//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// Graph InternalOutEdgeIterator - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec>
struct Iterator<Graph<Automaton<TAlphabet, TCargo, TGraphSpec> >, OutEdgeIterator>
{	
	typedef Iter<Graph<Automaton<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<OutEdgeIterator> > > Type;
};

template<typename TAlphabet, typename TCargo, typename TGraphSpec>
struct Iterator<Graph<Automaton<TAlphabet, TCargo, TGraphSpec> > const, OutEdgeIterator> :
    public Iterator<Graph<Automaton<TAlphabet, TCargo, TGraphSpec> >, OutEdgeIterator> {};


//////////////////////////////////////////////////////////////////////////////
// InternalOutEdgeIterator for Directed Graph
//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TSpec>
class Iter<Graph<Hmm<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > > 
{
public:
	typedef Graph<Hmm<TAlphabet, TCargo, TGraphSpec> > TGraph_;
	typedef typename EdgeDescriptor<TGraph_>::Type TEdgeDescriptor_;
	typedef typename VertexDescriptor<TGraph_>::Type TVertexDescriptor_;
	TGraph_ const* data_host;
	TVertexDescriptor_ data_source;
	TEdgeDescriptor_ data_edge;

	Iter()	
	{
		SEQAN_CHECKPOINT
	}
	
	Iter(TGraph_ const& _graph, TVertexDescriptor_ const v) : 
		data_host(&_graph),
		data_source(v),
		data_edge(getValue(_graph.data_model.data_vertex,v))
	{
		SEQAN_CHECKPOINT
	}
	
	Iter(Iter const& _iter) : 
		data_host(_iter.data_host),
		data_source(_iter.data_source),
		data_edge(_iter.data_edge)
	{
		SEQAN_CHECKPOINT
	}

	~Iter() {
		SEQAN_CHECKPOINT
	}

	Iter const&	operator = (Iter const & _other) {
		SEQAN_CHECKPOINT
		if (this == &_other) return *this;
		data_host = _other.data_host;
		data_source = _other.data_source;
		data_edge = _other.data_edge;
		return *this;
	}
//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
// Graph InternalOutEdgeIterator - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec>
struct Iterator<Graph<Hmm<TAlphabet, TCargo, TGraphSpec> >, OutEdgeIterator>
{	
	typedef Iter<Graph<Hmm<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<OutEdgeIterator> > > Type;
};

template<typename TAlphabet, typename TCargo, typename TGraphSpec>
struct Iterator<Graph<Hmm<TAlphabet, TCargo, TGraphSpec> > const, OutEdgeIterator> :
    public Iterator<Graph<Hmm<TAlphabet, TCargo, TGraphSpec> >, OutEdgeIterator> {};


//////////////////////////////////////////////////////////////////////////////
// Graph InternalOutEdgeIterator - General Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TIteratorSpec>
struct Value<Iter<TGraph, GraphIterator<InternalOutEdgeIterator<TIteratorSpec> > > >
{
	typedef typename EdgeDescriptor<TGraph>::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Value<Iter<TGraph const, GraphIterator<InternalOutEdgeIterator<TIteratorSpec> > > >
{
	typedef typename EdgeDescriptor<TGraph const>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TIteratorSpec>
struct Reference<Iter<TGraph, GraphIterator<InternalOutEdgeIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph, GraphIterator<InternalOutEdgeIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Reference<Iter<TGraph const, GraphIterator<InternalOutEdgeIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph const, GraphIterator<InternalOutEdgeIterator<TIteratorSpec> > > >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TIteratorSpec>
struct GetValue<Iter<TGraph, GraphIterator<InternalOutEdgeIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph, GraphIterator<InternalOutEdgeIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct GetValue<Iter<TGraph const, GraphIterator<InternalOutEdgeIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph const, GraphIterator<InternalOutEdgeIterator<TIteratorSpec> > > >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TIteratorSpec>
struct Spec<Iter<TGraph, GraphIterator<InternalOutEdgeIterator<TIteratorSpec> > > >
{
	typedef TIteratorSpec Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Spec<Iter<TGraph const, GraphIterator<InternalOutEdgeIterator<TIteratorSpec> > > >
{
	typedef TIteratorSpec Type;
};



//////////////////////////////////////////////////////////////////////////////
// Graph InternalOutEdgeIterator - Functions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline typename GetValue<Iter<Graph<Directed<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > > >::Type
getValue(Iter<Graph<Directed<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return it.data_edge;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TSpec>
inline typename GetValue<Iter<Graph<Hmm<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > > >::Type
getValue(Iter<Graph<Hmm<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return it.data_edge;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline typename GetValue<Iter<Graph<Tree<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > > >::Type
getValue(Iter<Graph<Tree<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return it.data_edge;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline typename GetValue<Iter<Graph<Undirected<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > > >::Type
getValue(Iter<Graph<Undirected<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return it.data_edge;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TSpec>
inline typename GetValue<Iter<Graph<Automaton<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > > >::Type
getValue(Iter<Graph<Automaton<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TGraphSpec> > TGraph;
	TGraph * g = const_cast<TGraph*>(it.data_host);
	return (&g->data_vertex[it.data_source].data_edge[it.data_pos]);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline typename Reference<Iter<Graph<Directed<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > > >::Type
value(Iter<Graph<Directed<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return it.data_edge;
}


//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TSpec>
inline typename Reference<Iter<Graph<Hmm<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > > >::Type
value(Iter<Graph<Hmm<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return it.data_edge;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline typename Reference<Iter<Graph<Tree<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > > >::Type
value(Iter<Graph<Tree<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return it.data_edge;
}


//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline typename Reference<Iter<Graph<Undirected<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > > >::Type
value(Iter<Graph<Undirected<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return it.data_edge;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TSpec>
inline typename Reference<Iter<Graph<Automaton<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > > >::Type
value(Iter<Graph<Automaton<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TGraphSpec> > TGraph;
	TGraph * g = const_cast<TGraph*>(it.data_host);
	return (&g->data_vertex[it.data_source].data_edge[it.data_pos]);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline typename Reference<Iter<TGraph, GraphIterator<InternalOutEdgeIterator<TSpec> > > >::Type
operator * (Iter<TGraph, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return value(it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline typename Host<Iter<TGraph, GraphIterator<InternalOutEdgeIterator<TSpec> > > >::Type const&
hostGraph(Iter<TGraph, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return *it.data_host;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline bool
atBegin(Iter<Graph<Directed<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (it.data_edge == getValue(it.data_host->data_vertex, it.data_source));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TSpec>
inline bool
atBegin(Iter<Graph<Hmm<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (it.data_edge == getValue(it.data_host->data_model.data_vertex, it.data_source));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline bool
atBegin(Iter<Graph<Tree<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (it.data_edge == getValue(it.data_host->data_vertex, it.data_source));
}


//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline bool
atBegin(Iter<Graph<Undirected<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (it.data_edge == getValue(it.data_host->data_vertex, it.data_source));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TSpec>
inline bool
atBegin(Iter<Graph<Automaton<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (it.data_begin == it.data_pos);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goBegin(Iter<Graph<Directed<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	it.data_edge = getValue(it.data_host->data_vertex,it.data_source);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goBegin(Iter<Graph<Hmm<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	it.data_edge = getValue(it.data_host->data_model.data_vertex,it.data_source);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goBegin(Iter<Graph<Tree<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	it.data_edge = getValue(it.data_host->data_vertex,it.data_source);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goBegin(Iter<Graph<Undirected<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	it.data_edge = getValue(it.data_host->data_vertex,it.data_source);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goBegin(Iter<Graph<Automaton<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	it.data_pos = it.data_begin;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline bool
atEnd(Iter<Graph<Directed<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (it.data_edge == 0);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TSpec>
inline bool
atEnd(Iter<Graph<Hmm<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (it.data_edge == 0);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline bool
atEnd(Iter<Graph<Tree<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (it.data_edge == 0);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline bool
atEnd(Iter<Graph<Undirected<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (it.data_edge == 0);
}

/////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TSpec>
inline bool
atEnd(Iter<Graph<Automaton<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (it.data_end == it.data_pos);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goEnd(Iter<Graph<Directed<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	it.data_edge = 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goEnd(Iter<Graph<Hmm<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	it.data_edge = 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goEnd(Iter<Graph<Tree<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	it.data_edge = 0;
}


//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goEnd(Iter<Graph<Undirected<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	it.data_edge = 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goEnd(Iter<Graph<Automaton<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	it.data_pos = it.data_end;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goNext(Iter<Graph<Directed<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	if (!atEnd(it)) it.data_edge = getNextT(it.data_edge);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goNext(Iter<Graph<Hmm<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	if (!atEnd(it)) it.data_edge = getNextT(it.data_edge);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goNext(Iter<Graph<Tree<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	if (!atEnd(it)) it.data_edge = getNextT(it.data_edge);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goNext(Iter<Graph<Undirected<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	if (!atEnd(it)) {
		if (it.data_source == getSource(it.data_edge)) it.data_edge = getNextS(it.data_edge);
		else it.data_edge = getNextT(it.data_edge);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goNext(Iter<Graph<Automaton<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TGraphSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	if (it.data_pos < it.data_end) ++it.data_pos;
	while (	(it.data_pos < it.data_end) &&
			(it.data_host->data_vertex[it.data_source].data_edge[it.data_pos].data_target == nilVal))
	{
				++it.data_pos;
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<InternalOutEdgeIterator<TSpec> > >&
operator ++(Iter<TGraph, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
SEQAN_CHECKPOINT
	goNext(it);
	return it;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<InternalOutEdgeIterator<TSpec> > >
operator ++(Iter<TGraph, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it, int)
{
	SEQAN_CHECKPOINT
	Iter<TGraph, GraphIterator<InternalOutEdgeIterator<TSpec> > > ret = it;
	goNext(it);
	return ret;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goPrevious(Iter<Graph<Directed<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	typedef typename EdgeType<Graph<Directed<TCargo, TGraphSpec> > >::Type TEdge;
	TEdge* current = getValue(it.data_host->data_vertex, it.data_source);
	if (current == it.data_edge) return;
	while ((current != 0) && (getNextT(current) != it.data_edge)) current = getNextT(current);
	it.data_edge = current;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goPrevious(Iter<Graph<Hmm<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	typedef typename EdgeType<Graph<Hmm<TAlphabet, TCargo, TGraphSpec> > >::Type TEdge;
	TEdge* current = getValue(it.data_host->data_model.data_vertex, it.data_source);
	if (current == it.data_edge) return;
	while ((current != 0) && (getNextT(current) != it.data_edge)) current = getNextT(current);
	it.data_edge = current;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goPrevious(Iter<Graph<Tree<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	typedef typename EdgeType<Graph<Tree<TCargo, TGraphSpec> > >::Type TEdge;
	TEdge* current = getValue(it.data_host->data_vertex, it.data_source);
	if (current == it.data_edge) return;
	while ((current != 0) && (getNextT(current) != it.data_edge)) current = getNextT(current);
	it.data_edge = current;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goPrevious(Iter<Graph<Undirected<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	typedef Graph<Undirected<TCargo, TGraphSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdge;
	TEdge* current = getValue(it.data_host->data_vertex, it.data_source);
	if (current == it.data_edge) return;
	while (current != 0) {
		if (it.data_source == getSource(current)) {
			if (it.data_edge == getNextS(current)) break;
			else current = getNextS(current);
		} else {
			if (it.data_edge == getNextT(current)) break;
			else current = getNextT(current);
		}
	}
	it.data_edge = current;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goPrevious(Iter<Graph<Automaton<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TGraphSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	if (it.data_pos != it.data_begin) --it.data_pos;
	while (	(it.data_pos != it.data_begin) &&
			(it.data_host->data_vertex[it.data_source].data_edge[it.data_pos].data_target == nilVal))
	{
				--it.data_pos;
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<InternalOutEdgeIterator<TSpec> > >&
operator --(Iter<TGraph, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	goPrevious(it);
	return it;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<InternalOutEdgeIterator<TSpec> > >
operator --(Iter<TGraph, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it, int)
{
	SEQAN_CHECKPOINT
	Iter<TGraph, GraphIterator<InternalOutEdgeIterator<TSpec> > > ret = it;
	goPrevious(it);
	return ret;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline bool
operator ==(Iter<Graph<Directed<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it1,
			Iter<Graph<Directed<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it2)
{
	SEQAN_CHECKPOINT
	return ((it1.data_edge==it2.data_edge) && 
			(it1.data_source==it2.data_source));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TSpec>
inline bool
operator ==(Iter<Graph<Hmm<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it1,
			Iter<Graph<Hmm<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it2)
{
	SEQAN_CHECKPOINT
	return ((it1.data_edge==it2.data_edge) && 
			(it1.data_source==it2.data_source));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline bool
operator ==(Iter<Graph<Tree<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it1,
			Iter<Graph<Tree<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it2)
{
	SEQAN_CHECKPOINT
	return ((it1.data_edge==it2.data_edge) && 
			(it1.data_source==it2.data_source));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline bool
operator ==(Iter<Graph<Undirected<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it1,
			Iter<Graph<Undirected<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it2)
{
	SEQAN_CHECKPOINT
	return ((it1.data_edge==it2.data_edge) && 
			(it1.data_source==it2.data_source));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TSpec>
inline bool
operator ==(Iter<Graph<Automaton<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it1,
			Iter<Graph<Automaton<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it2)
{
	SEQAN_CHECKPOINT
	return ((it1.data_pos==it2.data_pos) && 
			(it1.data_source==it2.data_source));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline bool
operator !=(Iter<Graph<Directed<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it1,
			Iter<Graph<Directed<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it2)
{
	SEQAN_CHECKPOINT
	return ((it1.data_edge!=it2.data_edge) || 
			(it1.data_source!=it2.data_source));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TSpec>
inline bool
operator !=(Iter<Graph<Hmm<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it1,
			Iter<Graph<Hmm<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it2)
{
	SEQAN_CHECKPOINT
	return ((it1.data_edge!=it2.data_edge) || 
			(it1.data_source!=it2.data_source));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline bool
operator !=(Iter<Graph<Tree<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it1,
			Iter<Graph<Tree<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it2)
{
	SEQAN_CHECKPOINT
	return ((it1.data_edge!=it2.data_edge) || 
			(it1.data_source!=it2.data_source));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline bool
operator !=(Iter<Graph<Undirected<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it1,
			Iter<Graph<Undirected<TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it2)
{
	SEQAN_CHECKPOINT
	return ((it1.data_edge!=it2.data_edge) || 
			(it1.data_source!=it2.data_source));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TSpec>
inline bool
operator !=(Iter<Graph<Automaton<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it1,
			Iter<Graph<Automaton<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it2)
{
	SEQAN_CHECKPOINT
	return ((it1.data_pos!=it2.data_pos) || 
			(it1.data_source!=it2.data_source));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline typename VertexDescriptor<TGraph>::Type 
sourceVertex(Iter<TGraph, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return it.data_source;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.Automaton#label
..class:Spec.Out-Edge Iterator
..class:Spec.Edge Iterator
..cat:Graph
..summary:Returns the label of the out-edge this iterator points to (for automatons).
..signature:label(it)
..param.it:An out-edge iterator.
...type:Spec.Out-Edge Iterator
...type:Spec.Edge Iterator
..returns:A label.
...type:Metafunction.Alphabet.
..remarks:The label function only works for out-edge iterators on automatons.
..include:seqan/graph_types.h
*/

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TSpec>
inline typename Alphabet<Graph<Automaton<TAlphabet, TCargo, TGraphSpec> > >::Type
label(Iter<Graph<Automaton<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	return TAlphabet(it.data_pos);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline typename VertexDescriptor<Graph<Directed<TCargo, TGraphSpec> > >::Type 
targetVertex(Iter<Graph<Directed<TCargo, TGraphSpec> > , GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return targetVertex(*it.data_host, it.data_edge);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TSpec>
inline typename VertexDescriptor<Graph<Hmm<TAlphabet, TCargo, TGraphSpec> > >::Type 
targetVertex(Iter<Graph<Hmm<TAlphabet, TCargo, TGraphSpec> > , GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return targetVertex(*it.data_host, it.data_edge);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline typename VertexDescriptor<Graph<Tree<TCargo, TGraphSpec> > >::Type 
targetVertex(Iter<Graph<Tree<TCargo, TGraphSpec> > , GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return targetVertex(*it.data_host, it.data_edge);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline typename VertexDescriptor<Graph<Undirected<TCargo, TGraphSpec> > >::Type 
targetVertex(Iter<Graph<Undirected<TCargo, TGraphSpec> > , GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	typedef Graph<Undirected<TCargo, TGraphSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TVertexDescriptor target = targetVertex(*it.data_host, it.data_edge);
	if (target != it.data_source) return target;
	else return sourceVertex(*it.data_host, it.data_edge);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TSpec>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TGraphSpec> > >::Type
targetVertex(Iter<Graph<Automaton<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return it.data_host->data_vertex[it.data_source].data_edge[it.data_pos].data_target;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
