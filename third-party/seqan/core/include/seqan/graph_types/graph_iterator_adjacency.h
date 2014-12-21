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

#ifndef SEQAN_HEADER_GRAPH_ITERATOR_ADJACENCY_H
#define SEQAN_HEADER_GRAPH_ITERATOR_ADJACENCY_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Graph AdjacencyIterator
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Adjacency Iterator:
..cat:Graph
..summary:Adjacency iterator for @Class.Graph@.
..signature:Iterator<TGraph, AdjacencyIterator>
..param.TGraph:A graph.
...type:Class.Graph
..general:Class.Iter
..see:Spec.Vertex Iterator
..see:Spec.Out-Edge Iterator
..see:Spec.Edge Iterator
..see:Spec.Bfs Iterator
..include:seqan/graph_types.h

.Memfunc.Adjacency Iterator#Iter
..summary:Constructor
..class:Spec.Adjacency Iterator
..signature:Iter()
..signature:Iter(graph, v)
..signature:Iter(iter)
..param.graph:The graph to iterate edges of.
...type:Class.Graph
..param.v:Vertex descriptor of vertex to iterate adjacent edges of.
...type:Metafunction.VertexDescriptor
..param.iter:Other iterator to copy-construct from.
...type:Spec.Adjacency Iterator
*/
template<typename TGraph, typename TSpec>
class Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > > 
{
public:
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor_;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator_;
	TOutEdgeIterator_ data_edge_it;

	Iter()	
	{
		SEQAN_CHECKPOINT
	}
	
	Iter(TGraph const& _graph, TVertexDescriptor_ const v) : 
		data_edge_it(_graph, v)
	{
		SEQAN_CHECKPOINT
	}
	
	~Iter() {
		SEQAN_CHECKPOINT
	}

	Iter(Iter const& _iter) : data_edge_it(_iter.data_edge_it)
	{
		SEQAN_CHECKPOINT
	}

	Iter const&	operator = (Iter const & _other) {
		SEQAN_CHECKPOINT
		if (this == &_other) return *this;
		data_edge_it = _other.data_edge_it;
		return *this;
	}
//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// Graph InternalAdjacencyIterator - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph>
struct Iterator<TGraph, AdjacencyIterator>
{	
	typedef Iter<TGraph, GraphIterator<InternalAdjacencyIterator<AdjacencyIterator> > > Type;
};

template<typename TGraph>
struct Iterator<TGraph const, AdjacencyIterator>
{	
	typedef Iter<TGraph const, GraphIterator<InternalAdjacencyIterator<AdjacencyIterator> > > Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TIteratorSpec>
struct Value<Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Value<Iter<TGraph const, GraphIterator<InternalAdjacencyIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph const, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TIteratorSpec>
struct Reference<Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TIteratorSpec> > > >
{
	typedef typename Reference<Iter<TGraph, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Reference<Iter<TGraph const, GraphIterator<InternalAdjacencyIterator<TIteratorSpec> > > >
{
	typedef typename Reference<Iter<TGraph const, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TIteratorSpec>
struct GetValue<Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TIteratorSpec> > > >
{
	typedef typename GetValue<Iter<TGraph, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct GetValue<Iter<TGraph const, GraphIterator<InternalAdjacencyIterator<TIteratorSpec> > > >
{
	typedef typename GetValue<Iter<TGraph const, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TIteratorSpec>
struct Spec<Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TIteratorSpec> > > >
{
	typedef TIteratorSpec Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Spec<Iter<TGraph const, GraphIterator<InternalAdjacencyIterator<TIteratorSpec> > > >
{
	typedef TIteratorSpec Type;
};

//////////////////////////////////////////////////////////////////////////////
// Graph InternalAdjacencyIterator - Functions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline typename GetValue<Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > > >::Type
getValue(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return targetVertex(it.data_edge_it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline typename GetValue<Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > > >::Type
value(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return getValue(it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline typename GetValue<Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > > >::Type
operator * (Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return value(it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline typename Host<Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > > >::Type const&
hostGraph(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return hostGraph(it.data_edge_it);
} 

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline bool
atBegin(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return atBegin(it.data_edge_it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline void
goBegin(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	goBegin(it.data_edge_it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline bool
atEnd(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (atEnd(it.data_edge_it));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline void
goEnd(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	goEnd(it.data_edge_it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline void
goNext(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	goNext(it.data_edge_it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >&
operator ++(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it)
{
SEQAN_CHECKPOINT
	goNext(it);
	return it;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >
operator ++(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it, int)
{
	SEQAN_CHECKPOINT
	Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > > ret = it;
	goNext(it);
	return ret;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline void
goPrevious(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	goPrevious(it.data_edge_it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >&
operator --(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	goPrevious(it);
	return it;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >
operator --(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it, int)
{
	SEQAN_CHECKPOINT
	Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > > ret = it;
	goPrevious(it);
	return ret;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline bool
operator ==(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it1,
			Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it2)
{
SEQAN_CHECKPOINT
	return (it1.data_edge_it==it2.data_edge_it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline bool
operator !=(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it1,
			Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it2)
{
SEQAN_CHECKPOINT
	return (it1.data_edge_it!=it2.data_edge_it);
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
