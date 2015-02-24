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

#ifndef SEQAN_INDEX_QGRAM_BUCKETREFINEMENT_H_
#define SEQAN_INDEX_QGRAM_BUCKETREFINEMENT_H_

//#define SEQAN_DEBUG

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TText>
class Index<TText, IndexSa<InfixSegment> >
{
public:
    typename Member<Index, FibreText>::Type       text;
    Holder<typename Fibre<Index, FibreSA>::Type>  sa;

    Index() {}

    Index(Index & other) :
        text(other.text),
        sa(other.sa)
    {}

    Index(Index const & other) :
        text(other.text),
        sa(other.sa)
    {}

    template <typename TText_>
    Index(TText_ & _text) :
        text(_text)
    {}

    template <typename TText_>
    Index(TText_ const & _text) :
        text(_text)
    {}
};

struct BucketRefinement_;
typedef Tag<BucketRefinement_> BucketRefinement;

template <typename TText, typename TShapeSpec>
class Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >:
    public Index<TText, IndexQGram<TShapeSpec> >
{
public:
    typedef Index<TText, IndexQGram<TShapeSpec> >    TBase;
    typedef Index<TText, IndexSa<InfixSegment> >     TIndexSa;

    TIndexSa    _indexSa;

    Index() :
        TBase()
    {
        _setHost(*this);
    }

    Index(Index & other) :
        TBase(static_cast<TBase &>(other)),
        _indexSa(other._indexSA)
    {
        _setHost(*this);
    }

    Index(Index const & other) :
        TBase(static_cast<TBase const &>(other)),
        _indexSa(other._indexSa)
    {
        _setHost(*this);
    }

    template <typename TText_>
    Index(TText_ & _text) :
        TBase(_text),
        _indexSa(_text)
    {
        _setHost(*this);
    }

    template <typename TText_>
    Index(TText_ const & _text) :
        TBase(_text),
        _indexSa(_text)
    {
        _setHost(*this);
    }

    template <typename TText_, typename TShape_>
    Index(TText_ & _text, TShape_ const & _shape) :
        TBase(_text, _shape),
        _indexSa(_text)
    {
        _setHost(*this);
    }

    template <typename TText_, typename TShape_>
    Index(TText_ const & _text, TShape_ const & _shape) :
        TBase(_text, _shape),
        _indexSa(_text)
    {
        _setHost(*this);
    }
};

// ============================================================================

template <typename TText, typename TShapeSpec, typename TSpec>
class Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > >
{
public:
    typedef Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >         TIndex;
    typedef Iter<TIndex, VSTree<TopDown<ParentLinks<TSpec> > > >            TParentLinksIter;

    typedef Index<TText, IndexQGram<TShapeSpec> >                           TTopIndex;
    typedef typename Iterator<TTopIndex, TopDown<TSpec> >::Type             TTopIterator;

    typedef Index<TText, IndexSa<InfixSegment> >                            TBottomIndex;
    typedef typename Iterator<TBottomIndex, TopDown<TSpec> >::Type          TBottomIterator;

    TTopIterator    _topIterator;
    TBottomIterator _bottomIterator;

    Iter() {}

    Iter(TIndex & _index) :
        _topIterator(),
        _bottomIterator(_index._indexSa)
    {
        _indexRequireTopDownIteration(_index);
        _topIterator = TTopIterator(static_cast<TTopIndex &>(_index));
        goRoot(_topIterator);
        goRoot(_bottomIterator);
    }

    Iter(Iter const & _origin) :
        _topIterator(_origin._topIterator),
        _bottomIterator(_origin._bottomIterator)
    {}

    Iter(TParentLinksIter const & _origin) :
        _topIterator(_origin._topIterator),
        _bottomIterator(_origin._bottomIterator)
    {}

    inline Iter const &
    operator=(Iter const & _origin)
    {
        _bottomIterator = _origin._bottomIterator;
        _topIterator = _origin._topIterator;
        return *this;
    }
};

// ============================================================================

template <typename TText, typename TShapeSpec, typename TSpec>
class Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<ParentLinks<TSpec> > > >
{
public:
    typedef Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >                 TIndex;
    typedef Iter<TIndex, VSTree<TopDown<TSpec> > >                                  TBase;

    typedef Index<TText, IndexQGram<TShapeSpec> >                                   TTopIndex;
    typedef typename Iterator<TTopIndex, TopDown<ParentLinks<TSpec> > >::Type       TTopIterator;

    typedef Index<TText, IndexSa<InfixSegment> >                                    TBottomIndex;
    typedef typename Iterator<TBottomIndex, TopDown<ParentLinks<TSpec> > >::Type    TBottomIterator;

    TTopIterator    _topIterator;
    TBottomIterator _bottomIterator;

    Iter() {}

    Iter(TIndex & _index) :
        _topIterator(),
        _bottomIterator(_index._indexSa)
    {
        _indexRequireTopDownIteration(_index);
        _topIterator = TTopIterator(static_cast<TTopIndex &>(_index));
        goRoot(_topIterator);
        goRoot(_bottomIterator);
    }

    Iter(Iter const & _origin) :
        _topIterator(_origin._topIterator),
        _bottomIterator(_origin._bottomIterator)
    {}

    inline Iter const &
    operator=(Iter const & _origin)
    {
        _bottomIterator = _origin._bottomIterator;
        _topIterator = _origin._topIterator;
        return *this;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TText, typename TShapeSpec>
struct Fibre<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreText>:
    public Fibre<Index<TText, IndexQGram<TShapeSpec> >, FibreText>
{};

template <typename TText, typename TShapeSpec>
struct Fibre<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreRawText>:
    public Fibre<Index<TText, IndexQGram<TShapeSpec> >, FibreRawText>
{};

template <typename TText, typename TShapeSpec>
struct Fibre<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreSA>:
    public Fibre<Index<TText, IndexQGram<TShapeSpec> >, FibreSA>
{};

template <typename TText, typename TShapeSpec>
struct Fibre<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreRawSA>:
    public Fibre<Index<TText, IndexQGram<TShapeSpec> >, FibreRawSA>
{};

template <typename TText, typename TShapeSpec>
struct Fibre<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreDir>:
    public Fibre<Index<TText, IndexQGram<TShapeSpec> >, FibreDir>
{};

template <typename TText, typename TShapeSpec>
struct Fibre<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreSADir>:
    public Fibre<Index<TText, IndexQGram<TShapeSpec> >, FibreSADir>
{};

template <typename TText, typename TShapeSpec>
struct Fibre<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreShape>:
    public Fibre<Index<TText, IndexQGram<TShapeSpec> >, FibreShape>
{};

template <typename TText, typename TShapeSpec>
struct Fibre<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreCounts>:
    public Fibre<Index<TText, IndexQGram<TShapeSpec> >, FibreCounts>
{};

template <typename TText, typename TShapeSpec>
struct Fibre<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreCountsDir>:
    public Fibre<Index<TText, IndexQGram<TShapeSpec> >, FibreCountsDir>
{};

template <typename TText, typename TShapeSpec>
struct Fibre<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreBucketMap>:
    public Fibre<Index<TText, IndexQGram<TShapeSpec> >, FibreBucketMap>
{};

// ============================================================================

template <typename TText, typename TShapeSpec>
struct DefaultIndexCreator<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreSA>:
    public DefaultIndexCreator<Index<TText, IndexSa<InfixSegment> >, FibreSA>
{};

// ============================================================================

template <typename TText, typename TShapeSpec, typename TSpec>
struct EdgeLabel<Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TSpec> > > :
    EdgeLabel<Iter<Index<TText, IndexQGram<TShapeSpec> >, VSTree<TSpec> > > {};

// ============================================================================

template <typename TText, typename TShapeSpec, typename TSpec>
struct Iterator<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > >
{
    typedef Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > >     Type;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TText>
inline typename Fibre<Index<TText, IndexSa<InfixSegment> >, FibreSA>::Type &
getFibre(Index<TText, IndexSa<InfixSegment> > & index, FibreSA)
{
    return value(index.sa);
}

template <typename TText>
inline typename Fibre<Index<TText, IndexSa<InfixSegment> > const, FibreSA>::Type &
getFibre(Index<TText, IndexSa<InfixSegment> > const & index, FibreSA)
{
    return value(index.sa);
}

// Works by creating the full SA and removing the suffixes of length < q.
template <typename TText, typename TShapeSpec>
inline bool indexCreate(Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > & index, FibreSADir, Default const)
{
    // Create QGram directory.
    resize(indexDir(index), _fullDirLength(index), Exact());
    createQGramIndexDirOnly(indexDir(index), indexBucketMap(index), indexText(index), indexShape(index), getStepSize(index));

    // Create full SA.
    indexCreate(index, FibreSA());

    // Remove too short suffixes from SA.
    _pruneSA(index);

    // Update indexSA host.
    _setHost(index);

    return true;
}

// Works by creating the q-gram directory and quick-sorting the SA bucket-wise.
//template <typename TText, typename TShapeSpec>
//inline bool indexCreate(Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > & index, FibreSADir, Default const)
//{
//    typedef Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >                 TIndex;
//    typedef Index<TText, IndexQGram<TShapeSpec> >                                   TTopIndex;
//
////    // Create QGram directory and refine suffix array.
//    indexCreate(static_cast<TTopIndex &>(index), FibreSADir(), Default());
//    _refineQGramIndex(indexSA(index), indexDir(index), indexText(index),
//                      weight(indexShape(index)), lengthSum(indexText(index)));
//
//    // Update indexSA host.
//    _setHost(index);
//
//    return true;
//}

template <typename TText, typename TShapeSpec>
void _setHost(Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > & index)
{
    setValue(index._indexSa.sa, indexSA(index));
}

template <typename TText, typename TShapeSpec>
void _pruneSA(Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > & index)
{
    typedef Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >     TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type                       TSA;
    typedef typename Iterator<TSA, Standard>::Type                      TSAIterator;

    TSA & sa = indexSA(index);

    TSAIterator saBegin = begin(sa, Standard());
    TSAIterator saEnd = end(sa, Standard());

    TSAIterator saOld = saBegin;
    TSAIterator saNew = saOld;

    while (saOld != saEnd)
    {
        if (suffixLength(*saOld, index) < weight(indexShape(index)))
        {
            ++saOld;
        }
        else
        {
            *saNew = *saOld;
            ++saOld;
            ++saNew;
        }
    }

    resize(sa, saNew - saBegin, Exact());
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline bool _implantSa(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > & it)
{
    //typedef Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >     TIndex;
    //typedef typename Value<TIndex>::Type                                TAlphabet;

    if (repLength(it._topIterator) < weight(indexShape(container(it._topIterator))))
        return false;

    SEQAN_ASSERT_EQ(repLength(it._bottomIterator), 0u);
    SEQAN_ASSERT(isRoot(it._bottomIterator));

    _historyPush(it._bottomIterator);

    value(it._bottomIterator).repLen = value(it._topIterator).repLen;
    value(it._bottomIterator).range = value(it._topIterator).range;
    value(it._bottomIterator).lastChar = value(it._topIterator).lastChar;

    return true;
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline bool _atTop(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > const & it)
{
    return isRoot(it._bottomIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline bool isRoot(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TSpec> > const & it)
{
    return isRoot(it._topIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline bool isLeaf(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TSpec> > const & it)
{
    // TODO(esiragusa): Fix isLeaf: last qgram is a top leaf but there is no subtree attached
    return isLeaf(it._topIterator) && isLeaf(it._bottomIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline void goRoot(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TSpec> > & it)
{
    goRoot(it._bottomIterator);
    goRoot(it._topIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline bool goDown(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > & it)
{
    //typedef Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >     TIndex;
    //typedef Pair<typename Size<TIndex>::Type>                           TSARange;

    if (_atTop(it))
    {
        if (goDown(it._topIterator))
            return true;

        if (!_implantSa(it))
            return false;
    }

    return goDown(it._bottomIterator);
}

// NOTE(esiragusa): I should have overloaded _goDownChar() instead.
template <typename TText, typename TShapeSpec, typename TSpec, typename TObject>
inline bool _goDownObject(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > & it,
                          TObject const & obj,
                          False const & /* tag */)
{
    if (_atTop(it))
    {
        if (_goDownObject(it._topIterator, obj, False()))
            return true;

        if (!_implantSa(it))
            return false;
    }

    return _goDownObject(it._bottomIterator, obj, False());
}

template <typename TText, typename TShapeSpec, typename TSpec, typename TString, typename TSize>
inline bool _goDownString(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > & it,
                          TString const & pattern,
                          TSize & lcp)
{
    TSize topLcp = 0;

    if (_atTop(it))
    {
        if (_goDownString(it._topIterator, pattern, lcp))
            return true;

        if (!_implantSa(it))
            return false;

        topLcp = lcp;
    }

    bool wentDown = _goDownString(it._bottomIterator, suffix(pattern, topLcp), lcp);

    lcp += topLcp;

    return wentDown;
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline bool goRight(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > & it)
{
    return _atTop(it) ? goRight(it._topIterator) : goRight(it._bottomIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline bool
goUp(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<ParentLinks<TSpec> > > > & it)
{
    if (repLength(it._bottomIterator) > repLength(it._topIterator))
    {
        SEQAN_ASSERT_NOT(isRoot(it._bottomIterator));

        goUp(it._bottomIterator);

        SEQAN_ASSERT_NOT(isRoot(it._bottomIterator));

        if (repLength(it._bottomIterator) <= repLength(it._topIterator))
            goRoot(it._bottomIterator);

        return true;
    }

    return goUp(it._topIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline typename Size<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > >::Type
countOccurrences(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TSpec> > const & it)
{
    return _atTop(it) ? countOccurrences(it._topIterator) : countOccurrences(it._bottomIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline typename SAValue<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > >::Type
getOccurrence(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TSpec> > const & it)
{
    return _atTop(it) ? getOccurrence(it._topIterator) : getOccurrence(it._bottomIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline typename Infix<typename Fibre<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreSA>::Type const>::Type
getOccurrences(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TSpec> > const & it)
{
    return _atTop(it) ? getOccurrences(it._topIterator) : getOccurrences(it._bottomIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline typename Size<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > >::Type
repLength(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > const & it)
{
    return _atTop(it) ? repLength(it._topIterator) : repLength(it._bottomIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline typename Infix<typename Fibre<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreText>::Type const>::Type
representative(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TSpec> > const & it)
{
    return _atTop(it) ? representative(it._topIterator) : representative(it._bottomIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline typename Size<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > >::Type
parentRepLength(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > const & it)
{
    return _atTop(it) ? parentRepLength(it._topIterator) : parentRepLength(it._bottomIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline typename Size<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > >::Type
parentEdgeLength(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > const & it)
{
    return _atTop(it) ? parentEdgeLength(it._topIterator) : parentEdgeLength(it._bottomIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline typename EdgeLabel<Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > >::Type
parentEdgeLabel(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > const & it)
{
    return _atTop(it) ? parentEdgeLabel(it._topIterator) : parentEdgeLabel(it._bottomIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline typename Value<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > >::Type
parentEdgeFirstChar(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > const & it)
{
    return _atTop(it) ? parentEdgeFirstChar(it._topIterator) : parentEdgeFirstChar(it._bottomIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline Pair<typename Size<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > >::Type>
range(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TSpec> > const & it)
{
    return _atTop(it) ? range(it._topIterator) : range(it._bottomIterator);
}

template <typename TText, typename TShapeSpec>
inline bool open(Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > & index, const char * fileName)
{
    typedef Index<TText, IndexQGram<TShapeSpec> >    TBaseIndex;

    if (open(static_cast<TBaseIndex &>(index), fileName))
    {
        setHost(index._indexSa, indexText(index));
        _setHost(index);
        return true;
    }

    return false;
}

template <typename TText, typename TShapeSpec>
inline bool save(Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > & index, const char * fileName)
{
    typedef Index<TText, IndexQGram<TShapeSpec> >    TBaseIndex;

    return save(static_cast<TBaseIndex &>(index), fileName);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INDEX_QGRAM_BUCKETREFINEMENT_H_
