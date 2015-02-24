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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// This file contains the Truncated VSTree iterator spec for IndexSa.
// ==========================================================================

#ifndef SEQAN_INDEX_SA_STREE_TRUNCATED_H_
#define SEQAN_INDEX_SA_STREE_TRUNCATED_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Truncated
// ----------------------------------------------------------------------------

template <typename TSpec = void>
struct Truncated;

// ----------------------------------------------------------------------------
// Class Iter
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TSpec>
class Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<Truncated<TSpec> > > > :
    public Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<TSpec> > >
{
public:
    typedef Index<TText, IndexSa<TIndexSpec> >          TIndex;
    typedef Iter<TIndex, VSTree<TopDown<TSpec> > >      TBase;
    typedef typename VertexDescriptor<TIndex>::Type     TVertexDesc;
    typedef typename Size<TIndex>::Type                 TDepth;

    TDepth depth;

    Iter() :
        TBase(),
        depth(MaxValue<TDepth>::VALUE)
    {}

    Iter(TIndex &_index) :
        TBase(_index),
        depth(MaxValue<TDepth>::VALUE)
    {}

    Iter(TIndex &_index, MinimalCtor) :
        TBase(_index, MinimalCtor()),
        depth(MaxValue<TDepth>::VALUE)
    {}

    Iter(TIndex &_index, TVertexDesc const &_vDesc) :
        TBase(_index, _vDesc),
        depth(MaxValue<TDepth>::VALUE)
    {}

    Iter(Iter const &_origin):
        TBase(_origin),
        depth(_origin.depth)
    {}

    inline Iter const &
    operator = (Iter const &_origin)
    {
        *static_cast<TBase*>(this) = _origin;
        depth = _origin.depth;

        return *this;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TText, typename TIndexSpec, typename TSpec, typename TDfsOrder>
inline bool _isLeaf(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<Truncated<TSpec> > > > const & it,
                    VSTreeIteratorTraits<TDfsOrder, True> const)
{
    typedef Index<TText, IndexSa<TIndexSpec> >                          TIndex;
    typedef Iter<TIndex, VSTree<TopDown<TSpec> > >                      TBase;

    return repLength(it) >= it.depth || _isLeaf(static_cast<TBase>(it), VSTreeIteratorTraits<TDfsOrder, True>());
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline bool isRightTerminal(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<Truncated<TSpec> > > > const & it)
{
    typedef Index<TText, IndexSa<TIndexSpec> >                          TIndex;
    typedef Iter<TIndex, VSTree<TopDown<TSpec> > >                      TBase;

    return repLength(it) >= it.depth || isRightTerminal(static_cast<TBase>(it));
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline typename Infix<typename Fibre<Index<TText, IndexSa<TIndexSpec> >, FibreSA>::Type const>::Type
getEmptyEdges(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<Truncated<TSpec> > > > const & it)
{
    typedef Index<TText, IndexSa<TIndexSpec> >                          TIndex;
    typedef Iter<TIndex, VSTree<TopDown<TSpec> > >                      TBase;

    if (repLength(it) >= it.depth)
        return getOccurrences(it);
    else
        return getEmptyEdges(static_cast<TBase>(it));
}

template < typename TIndex, class TSpec >
inline typename Infix< typename Fibre<TIndex, FibreSA>::Type const >::Type
getEmptyEdges(Iter< TIndex, VSTree<TSpec> > const &it)
{
    typedef typename Size<TIndex>::Type TSize;

    TIndex const &index = container(it);
    TSize repLen = repLength(it);
    Pair<TSize> saRange = range(it);

    TSize i2 = saRange.i1;
    while (i2 < saRange.i2 && suffixLength(saAt(i2, index), index) <= repLen)
        ++i2;

    return infix(indexSA(index), saRange.i1, i2);
}

}

#endif  // #ifndef SEQAN_INDEX_SA_STREE_TRUNCATED_H_
