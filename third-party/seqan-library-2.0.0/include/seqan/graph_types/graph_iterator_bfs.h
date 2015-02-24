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

#ifndef SEQAN_HEADER_GRAPH_ITERATOR_BFS_H
#define SEQAN_HEADER_GRAPH_ITERATOR_BFS_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph BfsIterator
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class BfsIterator
 * @extends Iter
 * @headerfile <seqan/graph_types.h>
 * @brief Iterate vertices of a graph in breadth-first fashion.
 *
 * @signature Iterator<TGraph, BfsIterator>::Type;
 * @signature template <typename TGraph, typename TSpec>
 *            class Iter<TGraph, GraphIterator<InternalBfsIterator<TSpec> >;
 *
 * @tparam TGraph The graph to iterate the vertices of.
 *
 * The first signature is the signature of the corresponding @link ContainerConcept#Iterator graph's Iterator @endlink
 * metafunction call.  The second call is the internal definition of the type.  You should always get this type using
 * the metafunction call from the first signature.
 *
 *
 * @fn BfsIterator::BfsIterator
 * @brief Constructor.
 *
 * @signature Iter::Iter();
 * @signature Iter::Iter(iter);
 * @signature Iter::Iter(graph, v);
 *
 * @param[in] iter  Other BfsIterator to copy from.
 * @param[in] graph The @link Graph @endlink to iterate vertices in BFS fashion.
 * @param[in] v     The descriptor of the vertex to start BFS iteration.
 */

template<typename TGraph, typename TSpec>
class Iter<TGraph, GraphIterator<InternalBfsIterator<TSpec> > >
{
public:
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor_;
    TGraph const* data_host;
    TVertexDescriptor_ data_source;
    String<bool> data_tokenMap;
    std::deque<TVertexDescriptor_> data_queue;

    void _init() {
        resizeVertexMap(data_tokenMap, *data_host);
        typedef typename Iterator<String<bool>, Rooted>::Type TIter;
        TIter it = begin(data_tokenMap);
        for(;!atEnd(it);goNext(it)) {
            assignValue(it,false);
        }
        assignProperty(data_tokenMap, data_source, true);
        data_queue.clear();
        data_queue.push_back(data_source);
    }

    Iter()
    {
        SEQAN_CHECKPOINT
    }

    Iter(TGraph const& _graph, TVertexDescriptor_ v) :
        data_host(&_graph),
        data_source(v)
    {
        SEQAN_CHECKPOINT
        _init();
    }


    ~Iter() {
        SEQAN_CHECKPOINT
    }

    Iter(Iter const& _iter) :
        data_host(_iter.data_host),
        data_source(_iter.data_source),
        data_tokenMap(_iter.data_tokenMap),
        data_queue(_iter.data_queue)
    {
        SEQAN_CHECKPOINT
    }

    Iter const&    operator = (Iter const & _other) {
        SEQAN_CHECKPOINT
        if (this == &_other) return *this;
        data_host=_other.data_host;
        data_source=_other.data_source;
        data_tokenMap=_other.data_tokenMap;
        data_queue=_other.data_queue;
        return *this;
    }
//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// Graph InternalBfsIterator - Metafunctions
//////////////////////////////////////////////////////////////////////////////
template<typename TGraph>
struct Iterator<TGraph, BfsIterator>
{
    typedef Iter<TGraph, GraphIterator<InternalBfsIterator<BfsIterator> > > Type;
};

template<typename TGraph>
struct Iterator<TGraph const, BfsIterator>
{
    typedef Iter<TGraph const, GraphIterator<InternalBfsIterator<BfsIterator> > > Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Value<Iter<TGraph, GraphIterator<InternalBfsIterator<TIteratorSpec> > > >
{
    typedef typename Value<Iter<TGraph, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Value<Iter<TGraph const, GraphIterator<InternalBfsIterator<TIteratorSpec> > > >
{
    typedef typename Value<Iter<TGraph const, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Reference<Iter<TGraph, GraphIterator<InternalBfsIterator<TIteratorSpec> > > >
{
    typedef typename Reference<Iter<TGraph, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Reference<Iter<TGraph const, GraphIterator<InternalBfsIterator<TIteratorSpec> > > >
{
    typedef typename Reference<Iter<TGraph const, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct GetValue<Iter<TGraph, GraphIterator<InternalBfsIterator<TIteratorSpec> > > >
{
    typedef typename GetValue<Iter<TGraph, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct GetValue<Iter<TGraph const, GraphIterator<InternalBfsIterator<TIteratorSpec> > > >
{
    typedef typename GetValue<Iter<TGraph const, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////
// Graph InternalBfsIterator - Functions
//////////////////////////////////////////////////////////////////////////////
template<typename TGraph, typename TSpec>
inline typename GetValue<Iter<TGraph, GraphIterator<InternalBfsIterator<TSpec> > > >::Type
getValue(Iter<TGraph, GraphIterator<InternalBfsIterator<TSpec> > >& it)
{
    SEQAN_CHECKPOINT
    return it.data_queue.front();
}

template<typename TGraph, typename TSpec>
inline typename GetValue<Iter<TGraph, GraphIterator<InternalBfsIterator<TSpec> > > >::Type
value(Iter<TGraph, GraphIterator<InternalBfsIterator<TSpec> > >& it)
{
    SEQAN_CHECKPOINT
    return getValue(it);
}

template<typename TGraph, typename TSpec>
inline typename GetValue<Iter<TGraph, GraphIterator<InternalBfsIterator<TSpec> > > >::Type
operator * (Iter<TGraph, GraphIterator<InternalBfsIterator<TSpec> > >& it)
{
    SEQAN_CHECKPOINT
    return value(it);
}

template<typename TGraph, typename TSpec>
inline typename Host<Iter<TGraph, GraphIterator<InternalBfsIterator<TSpec> > > >::Type const&
hostGraph(Iter<TGraph, GraphIterator<InternalBfsIterator<TSpec> > >& it)
{
    SEQAN_CHECKPOINT
    return *it.data_host;
}

template<typename TGraph, typename TSpec>
inline bool
atBegin(Iter<TGraph, GraphIterator<InternalBfsIterator<TSpec> > >& it)
{
    SEQAN_CHECKPOINT
    if (it.data_queue.empty()) return false;
    else return (it.data_queue.front() == it.data_source);
}

template<typename TGraph, typename TSpec>
inline void
goBegin(Iter<TGraph, GraphIterator<InternalBfsIterator<TSpec> > >& it)
{
    SEQAN_CHECKPOINT
    it._init();
}

template<typename TGraph, typename TSpec>
inline bool
atEnd(Iter<TGraph, GraphIterator<InternalBfsIterator<TSpec> > >& it)
{
    SEQAN_CHECKPOINT
    return (it.data_queue.empty());
}

template<typename TGraph, typename TSpec>
inline void
goEnd(Iter<TGraph, GraphIterator<InternalBfsIterator<TSpec> > >& it)
{
    SEQAN_CHECKPOINT
    it.data_queue.clear();
}

template<typename TGraph, typename TSpec>
inline void
goNext(Iter<TGraph, GraphIterator<InternalBfsIterator<TSpec> > >& it)
{
    SEQAN_CHECKPOINT
    if (it.data_queue.empty()) return;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    TVertexDescriptor u = it.data_queue.front();
    it.data_queue.pop_front();
    typedef typename Iterator<TGraph, AdjacencyIterator>::Type TAdjacencyIterator;
    TAdjacencyIterator itad(*it.data_host,u);
    for(;!atEnd(itad);goNext(itad)) {
        TVertexDescriptor v = getValue(itad);
        if (getProperty(it.data_tokenMap, v) == false) {
            assignProperty(it.data_tokenMap, v, true);
            it.data_queue.push_back(v);
        }
    }
}

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<InternalBfsIterator<TSpec> > >&
operator ++(Iter<TGraph, GraphIterator<InternalBfsIterator<TSpec> > >& it)
{
SEQAN_CHECKPOINT
    goNext(it);
    return it;
}

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<InternalBfsIterator<TSpec> > >
operator ++(Iter<TGraph, GraphIterator<InternalBfsIterator<TSpec> > >& it, int)
{
    SEQAN_CHECKPOINT
    Iter<TGraph, GraphIterator<InternalBfsIterator<TSpec> > > ret = it;
    goNext(it);
    return ret;
}

template<typename TGraph, typename TSpec>
inline bool
operator ==(Iter<TGraph, GraphIterator<InternalBfsIterator<TSpec> > >& it1,
            Iter<TGraph, GraphIterator<InternalBfsIterator<TSpec> > >& it2)
{
    SEQAN_CHECKPOINT
    return ((it1.data_source==it2.data_source) &&
            (it1.data_tokenMap==it2.data_tokenMap) &&
            (it1.data_queue==it2.data_queue));
}

template<typename TGraph, typename TSpec>
inline bool
operator !=(Iter<TGraph, GraphIterator<InternalBfsIterator<TSpec> > >& it1,
            Iter<TGraph, GraphIterator<InternalBfsIterator<TSpec> > >& it2)
{
    SEQAN_CHECKPOINT
    return ((it1.data_source!=it2.data_source) ||
            (it1.data_tokenMap!=it2.data_tokenMap) ||
            (it1.data_queue!=it2.data_queue));
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
