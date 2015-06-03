// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
// Copyright (c) 2013 NVIDIA Corporation
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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

//SEQAN_NO_DDDOC:do not generate documentation for this file

#ifndef INDEX_FM_STREE_H_
#define INDEX_FM_STREE_H_

namespace seqan {

// ==========================================================================
// Classes
// ==========================================================================

// ----------------------------------------------------------------------------
// Class HistoryStackFM_
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
struct HistoryStackFM_
{
    Pair<TSize> range;
    TSize       repLen;
    TAlphabet   lastChar;

    SEQAN_HOST_DEVICE
    HistoryStackFM_() {}

    template <typename TSize_, typename TAlphabet_>
    SEQAN_HOST_DEVICE
    HistoryStackFM_(Pair<TSize_> const &_range, TSize_ _repLen, TAlphabet_ _lastChar):
        range(_range),
        repLen(_repLen),
        lastChar(_lastChar)
    {}

    SEQAN_HOST_DEVICE
    HistoryStackFM_ const &
    operator=(HistoryStackFM_ const & _origin)
    {
        range = _origin.range;
        repLen = _origin.repLen;
        lastChar = _origin.lastChar;
        return _origin;
    }
};

// ----------------------------------------------------------------------------
// Class VertexFM
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
struct VertexFM
{
    Pair<TSize> range;
    TSize       repLen;
    TAlphabet   lastChar;

    SEQAN_HOST_DEVICE
    VertexFM() :
        range(0, 0),
        repLen(0),
        lastChar(0)
    {}

    SEQAN_HOST_DEVICE
    VertexFM(MinimalCtor) :
        range(0, 0),
        repLen(0),
        lastChar(0)
    {}

    SEQAN_HOST_DEVICE
    VertexFM(Pair<TSize> newCurrentRange, TSize newRepLen, TAlphabet newChar) :
        range(newCurrentRange),
        repLen(newRepLen),
        lastChar(newChar)
    {}

    SEQAN_HOST_DEVICE
    VertexFM(VertexFM const & other) :
        range(other.range),
        repLen(other.repLen),
        lastChar(other.lastChar)
    {}

    SEQAN_HOST_DEVICE inline
    VertexFM &
    operator = (VertexFM const & _origin)
    {
        range = _origin.range;
        repLen = _origin.repLen;
        lastChar = _origin.lastChar;
        return *this;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction VertexDescriptor                                        [Index]
// ----------------------------------------------------------------------------

template < typename TText, typename TOccSpec, typename TIndexSpec>
struct VertexDescriptor<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >
{
    typedef Index<TText,FMIndex<TOccSpec, TIndexSpec> >     TIndex_;
    typedef typename Size<TIndex_>::Type                    TSize_;
    typedef typename Value<TIndex_>::Type                   TAlphabet_;

    typedef VertexFM<TSize_, TAlphabet_>                    Type;
};

// ----------------------------------------------------------------------------
// Metafunction HistoryStackEntry_                                      [Index]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TSpec, typename TIterSpec>
struct HistoryStackEntry_<Iter<Index<TText, FMIndex<TOccSpec, TSpec> >,
                               VSTree< TopDown< ParentLinks<TIterSpec> > > > >
{
    typedef HistoryStackFM_<typename Size<Index<TText, FMIndex<TOccSpec, TSpec> > >::Type,
                             typename Value<Index<TText, FMIndex<TOccSpec, TSpec> > >::Type>    Type;
};

// ----------------------------------------------------------------------------
// Metafunction EdgeLabel                                            [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
struct EdgeLabel<Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TSpec> > >
{
    typedef typename Value<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >::Type Type;
};


// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _indexRequireTopDownIteration()                             [Index]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TIndexSpec>
SEQAN_HOST_DEVICE inline void _indexRequireTopDownIteration(Index<TText, FMIndex<TOccSpec, TIndexSpec> > & index)
{
    indexRequire(index, FibreSALF());
}

// ----------------------------------------------------------------------------
// Function begin()                                                  [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
inline
typename Iterator<Index<TText,FMIndex<TOccSpec, TIndexSpec> >, TSpec>::Type
begin(Index<TText, FMIndex<TOccSpec, TIndexSpec> > & index, TSpec const /*Tag*/)
{
    typedef typename Iterator<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, TSpec>::Type TIter;

    return TIter(index);
}

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
inline
typename Iterator<Index<TText,FMIndex<TOccSpec, TIndexSpec> > const, TSpec>::Type
begin(Index<TText, FMIndex<TOccSpec, TIndexSpec> > const & index, TSpec const /*Tag*/)
{
    typedef typename Iterator<Index<TText, FMIndex<TOccSpec, TIndexSpec> > const, TSpec>::Type TIter;

    return TIter(index);
}

// ----------------------------------------------------------------------------
// Function _isRoot()                                                [Iterator]
// ----------------------------------------------------------------------------

template <typename TAlphabet, typename TSize>
SEQAN_HOST_DEVICE inline bool _isRoot(VertexFM<TAlphabet, TSize> const & value)
{
    return _isSizeInval(value.range.i2);
}

// ----------------------------------------------------------------------------
// Function _isLeaf()                                                [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder>
SEQAN_HOST_DEVICE inline bool
_isLeaf(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TSpec> > const & it,
        VSTreeIteratorTraits<TDfsOrder, False> const)
{
    return value(it).range.i1 >= value(it).range.i2;
}

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder>
SEQAN_HOST_DEVICE inline bool
_isLeaf(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TSpec> > const & it,
        VSTreeIteratorTraits<TDfsOrder, True> const)
{
    return _isLeaf(it, VSTreeIteratorTraits<TDfsOrder, False>()) &&
           isSentinel(indexLF(container(it)), value(it).range.i1);
}

// ----------------------------------------------------------------------------
// Function _getNodeByChar()                                         [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TChar>
SEQAN_HOST_DEVICE inline bool
_getNodeByChar(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it,
               typename VertexDescriptor<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >::Type const & vDesc,
               Pair<typename Size<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >::Type> & _range,
               TChar c)
{
    typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> >        TIndex;
    typedef typename Fibre<TIndex, FibreLF>::Type               TLF;

    TIndex const & index = container(it);
    TLF const & lf = indexLF(index);

    _range = range(index, vDesc);

    _range.i1 = lf(_range.i1, c);
    _range.i2 = lf(_range.i2, c);

    return _range.i1 < _range.i2;
}

// ----------------------------------------------------------------------------
// Function _goDownChar()                                            [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TChar>
SEQAN_HOST_DEVICE inline bool
_goDownChar(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > & it, TChar c)
{
    typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> >        TIndex;
    typedef Pair<typename Size<TIndex>::Type>                   TRange;

    TRange _range;

    if (_getNodeByChar(it, value(it), _range, c))
    {
        _historyPush(it);

        value(it).range = _range;
        value(it).lastChar = c;
        value(it).repLen++;

        return true;
    }

    return false;
}

// ----------------------------------------------------------------------------
// Function _goDown()                                                [Iterator]
// ----------------------------------------------------------------------------

// TODO(esiragusa): Implement this.
template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder>
SEQAN_HOST_DEVICE inline bool
_goDown(Iter<Index<TText, FMIndex<TOccSpec,TIndexSpec> >, VSTree<TopDown<TSpec> > > & it,
        VSTreeIteratorTraits<TDfsOrder, False> const)
{
    return _goDown(it, VSTreeIteratorTraits<TDfsOrder, True>());
}

template <typename TText, typename TOccSpec, typename TSpec, typename TIndexSpec, typename TDfsOrder>
SEQAN_HOST_DEVICE inline bool
_goDown(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > & it,
        VSTreeIteratorTraits<TDfsOrder, True> const)
{
    typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> >    TIndex;
    typedef typename Value<TIndex>::Type                    TAlphabet;
//    typedef typename ValueSize<TAlphabet>::Type             TAlphabetSize;

    // NOTE(esiragusa): isLeaf() early exit is slower on CUDA.
    // NOTE(esiragusa): this should be faster only for texts over small alphabets consisting of few/long sequences.
#ifndef __CUDA_ARCH__
    if (isLeaf(it)) return false;
#endif

    // TODO(esiragusa): Fix increment for alphabets with qualities.
//    for (TAlphabetSize c = 0; c < ValueSize<TAlphabet>::VALUE; ++c)
    for (TAlphabet c = 0; ordValue(c) < ValueSize<TAlphabet>::VALUE; ++c)
        if (_goDownChar(it, c)) return true;

    return false;
}

// ----------------------------------------------------------------------------
// Function _goDownString()                                          [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TString, typename TSize>
SEQAN_HOST_DEVICE inline bool
_goDownString(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > &it,
              TString const & string,
              TSize & lcp)
{
    typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> >        TIndex;
    typedef Pair<typename Size<TIndex>::Type>                   TRange;
    typedef typename Iterator<TString const, Standard>::Type    TStringIter;

    _historyPush(it);

    TStringIter stringIt = begin(string, Standard());
    TStringIter stringEnd = end(string, Standard());

    for (lcp = 0; stringIt != stringEnd; ++stringIt, ++lcp)
    {
        TRange _range;

        // NOTE(esiragusa): isLeaf() early exit is slower on CUDA.
        // NOTE(esiragusa): this should be faster only for texts over small alphabets consisting of few/long sequences.
#ifdef __CUDA_ARCH__
        if (!_getNodeByChar(it, value(it), _range, value(stringIt))) break;
#else
        if (isLeaf(it) || !_getNodeByChar(it, value(it), _range, value(stringIt))) break;
#endif

        value(it).range = _range;
    }

    value(it).repLen += lcp;

    if (lcp) value(it).lastChar = value(stringIt - 1);

    return stringIt == stringEnd;
}

// ----------------------------------------------------------------------------
// Function _goRight()                                               [Iterator]
// ----------------------------------------------------------------------------

// TODO(esiragusa): Implement this.
template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder>
SEQAN_HOST_DEVICE inline bool
_goRight(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > & it,
         VSTreeIteratorTraits<TDfsOrder, False> const)
{
    return _goRight(it, VSTreeIteratorTraits<TDfsOrder, True>());
}

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder>
SEQAN_HOST_DEVICE inline bool
_goRight(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > & it,
         VSTreeIteratorTraits<TDfsOrder, True> const)
{
    typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> >        TIndex;
    typedef typename Value<TIndex>::Type                        TAlphabet;
    typedef Pair<typename Size<TIndex>::Type>                   TRange;

    typedef typename VertexDescriptor<TIndex>::Type             TVertexDescriptor;

    if (isRoot(it)) return false;

    TVertexDescriptor parentDesc = nodeUp(it);
    TRange _range;

    // TODO(esiragusa): Fix increment for alphabets with qualities.
//    for (TAlphabetSize c = ordValue(value(it).lastChar) + 1; c < ValueSize<TAlphabet>::VALUE; ++c)
    for (value(it).lastChar++; ordValue(value(it).lastChar) < ValueSize<TAlphabet>::VALUE; value(it).lastChar++)
    {
        if (_getNodeByChar(it, parentDesc, _range, value(it).lastChar))
        {
            value(it).range = _range;

            return true;
        }
    }

    return false;
}

// ----------------------------------------------------------------------------
// Function _goUp()                                                  [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
SEQAN_HOST_DEVICE inline bool
_goUp(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > & it)
{
    if (isRoot(it)) return false;

    _historyPop(it);
    return true;
}

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
SEQAN_HOST_DEVICE inline bool
_goUp(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > & it)
{
    if (isRoot(it)) return false;

    _historyPop(it);
    return true;
}

// ----------------------------------------------------------------------------
// Function nodeUp()                                                 [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TIndexSpec, class TSpec>
SEQAN_HOST_DEVICE inline typename VertexDescriptor<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >::Type
nodeUp(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree< TopDown< ParentLinks<TSpec> > > > const & it)
{
    typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> >    TIndex;
    typedef typename VertexDescriptor<TIndex>::Type         TVertexDescriptor;

    if (!empty(it.history))
        return TVertexDescriptor(back(it.history).range, back(it.history).repLen, back(it.history).lastChar);
    else
        return value(it);
}

// ----------------------------------------------------------------------------
// Function _historyPush()                                           [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
SEQAN_HOST_DEVICE inline void
_historyPush(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec > >, VSTree<TopDown<TSpec> > > & it)
{
    it._parentDesc = value(it);
}

template < typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
SEQAN_HOST_DEVICE inline void
_historyPush(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec > >, VSTree<TopDown<ParentLinks<TSpec> > > > & it)
{
    typedef Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec > >, VSTree<TopDown<ParentLinks<TSpec> > > > TIter;

    typename HistoryStackEntry_<TIter>::Type h;

    h.range = value(it).range;
    h.repLen = value(it).repLen;
    h.lastChar = value(it).lastChar;

    appendValue(it.history, h);
}

// ----------------------------------------------------------------------------
// Function _historyPop()                                            [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
SEQAN_HOST_DEVICE inline void
_historyPop(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec > >, VSTree<TopDown<TSpec> > > & it)
{
    value(it).range = back(it.history).range;
    value(it).repLen = back(it.history).repLen;
    value(it).lastChar = back(it.history).lastChar;
}

template < typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
SEQAN_HOST_DEVICE inline void
_historyPop(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec > >, VSTree<TopDown<ParentLinks<TSpec> > > > & it)
{
    value(it).range = back(it.history).range;
    value(it).repLen = back(it.history).repLen;
    value(it).lastChar = back(it.history).lastChar;
    eraseBack(it.history);
}

// ----------------------------------------------------------------------------
// Function repLength()                                              [Iterator]
// ----------------------------------------------------------------------------

template <typename TIndex, typename TAlphabet, typename TSize>
SEQAN_HOST_DEVICE inline typename Size<TIndex>::Type
repLength(TIndex const &, VertexFM<TAlphabet, TSize> const & vDesc)
{
    return vDesc.repLen;
}

// ----------------------------------------------------------------------------
// Function parentEdgeLabel()                                        [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
SEQAN_HOST_DEVICE inline typename EdgeLabel<Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec > >, VSTree<TopDown<TSpec> > > >::Type
parentEdgeLabel(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec > >, VSTree<TopDown<TSpec> > > const & it)
{
    return value(it).lastChar;
}

// ----------------------------------------------------------------------------
// Function parentEdgeFirstChar()                                    [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
SEQAN_HOST_DEVICE inline typename Value<Index<TText, FMIndex<TOccSpec, TIndexSpec > > >::Type
parentEdgeFirstChar(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec > >, VSTree<TopDown<TSpec> > > const & it)
{
    return value(it).lastChar;
}

}
#endif  // INDEX_FM_STREE_H_
