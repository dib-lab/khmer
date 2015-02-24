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
// TopDown Iterators for the QGram Index.
// ==========================================================================
// Define SEQAN_INDEX_QGRAM_TREE to let the edge labels be tree-like infixes
// of arbitrary length instead of trie-like unitary symbols.
// ==========================================================================

#ifndef SEQAN_INDEX_QGRAM_STREE_H_
#define SEQAN_INDEX_QGRAM_STREE_H_

//#define SEQAN_DEBUG

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TSize, typename TShape>
struct VertexQGram
{
    typedef typename Host<TShape>::Type                 TAlphabet;
    typedef typename Value<TShape>::Type                THashValue;

    Pair<TSize>         range;
    Pair<THashValue>    hash;
    TSize               repLen;
    TAlphabet           lastChar;

    VertexQGram() :
        range(0, 0),
        hash(0, 0),
        repLen(0),
        lastChar(0)
    {}

    VertexQGram(MinimalCtor) :
        range(0, 0),
        hash(0, 0),
        repLen(0),
        lastChar(0)
    {}

    VertexQGram(VertexQGram const & other) :
        range(other.range),
        hash(other.hash),
        repLen(other.repLen),
        lastChar(other.lastChar)
    {}
};

template <typename TSize, typename TShape>
struct HistoryStackQGram_
{
    typedef typename Host<TShape>::Type                 TAlphabet;
    typedef typename Value<TShape>::Type                THashValue;

    Pair<TSize>         range;
    Pair<THashValue>    hash;
    TAlphabet           lastChar;

    HistoryStackQGram_() {}
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TText, typename TShapeSpec, typename TIndexSpec>
struct VertexDescriptor<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> > >
{
private:
    typedef Index<TText, IndexQGram<TIndexSpec> >       TIndex;
    typedef typename Fibre<TIndex, QGramShape>::Type    TShape;
    typedef typename Size<TIndex>::Type                 TSize;

public:
    typedef VertexQGram<TSize, TShape>                  Type;
};

template <typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec>
struct HistoryStackEntry_<Iter<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >,
                               VSTree<TopDown<ParentLinks<TSpec> > > > >
{
private:
    typedef Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >   TIndex;
    typedef typename Fibre<TIndex, QGramShape>::Type            TShape;
    typedef typename Size<TIndex>::Type                         TSize;

public:
    typedef HistoryStackQGram_<TSize, TShape>                   Type;
};

#ifndef SEQAN_INDEX_QGRAM_TREE
template <typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec>
struct EdgeLabel<Iter<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >, VSTree<TSpec> > >
{
    typedef typename Value<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> > >::Type Type;
};
#endif

template <typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec>
struct Iterator<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >, VSTree<TSpec> >
{
    typedef Iter<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >, VSTree<TSpec> >    Type;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TText, typename TShapeSpec, typename TIndexSpec>
void _indexRequireTopDownIteration(Index<TText, IndexQGram<TShapeSpec, TIndexSpec> > & index)
{
    indexRequire(index, QGramSADir());
}

template <typename TSize, typename TShape>
inline bool _isRoot(VertexQGram<TSize, TShape> const & value)
{
    return _isSizeInval(value.range.i2);
}

template <typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec>
inline bool isLeaf(Iter<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >, VSTree<TSpec> > const & it)
{
    return value(it).hash.i1 + 1 >= value(it).hash.i2;
}

template <typename TIndex, typename TSize, typename TShape>
inline typename Size<TIndex>::Type
repLength(TIndex const &, VertexQGram<TSize, TShape> const & vDesc)
{
    return vDesc.repLen;
}

#ifndef SEQAN_INDEX_QGRAM_TREE
template <typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec>
inline typename EdgeLabel<Iter<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > >::Type
parentEdgeLabel(Iter<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it)
{
    return value(it).lastChar;
}

template <typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec>
inline typename Value<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> > >::Type
parentEdgeFirstChar(Iter<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it)
{
    return value(it).lastChar;
}
#endif

template <typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec>
inline void goRoot(Iter<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >, VSTree<TSpec> > & it)
{
    typedef Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >   TIndex;
    typedef typename Fibre<TIndex, QGramShape>::Type            TShape;

    _historyClear(it);
    clear(it);
    if (!empty(indexSA(container(it))))
        _setSizeInval(value(it).range.i2);

    value(it).hash.i1 = 0;
    value(it).hash.i2 = ValueSize<TShape>::VALUE;
    value(it).repLen = 0;

#ifdef SEQAN_DEBUG
    std::cout << "[" << value(it).hash.i1 << "," << value(it).hash.i2 << ")" << std::endl;
#endif
}

template <typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec>
inline void
_historyPush(Iter<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > & it)
{
    it._parentDesc = value(it);
}

template <typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec>
inline void
_historyPush(Iter<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > & it)
{
    typedef Iter<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > TIter;
    typename HistoryStackEntry_<TIter>::Type h;
    h.range = value(it).range;
    h.hash = value(it).hash;
    h.lastChar = value(it).lastChar;

    appendValue(it.history, h);
}

template <typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec, typename TValue>
inline bool
_getNodeByChar(Iter<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >, VSTree<TSpec> > const & it,
               TValue c,
               typename VertexDescriptor<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> > >::Type & childDesc)
{
    typedef Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >   TIndex;
    typedef typename Fibre<TIndex, QGramDir>::Type              TDir;
    typedef typename Fibre<TIndex, QGramShape>::Type            TShape;
    typedef typename Host<TShape>::Type                         TAlphabet;
    typedef typename Value<TShape>::Type                        THValue;

    if (isLeaf(it))
        return false;

    TIndex const & index = container(it);
    TDir const & dir = indexDir(index);
    TShape const & shape = indexShape(index);

    childDesc = value(it);
    childDesc.lastChar = c;

    // TODO(esiragusa): Remove call to pow
    THValue h = static_cast<THValue>(pow(static_cast<double>(ValueSize<TAlphabet>::VALUE), static_cast<int>(weight(shape) - childDesc.repLen - 1)));
    childDesc.hash.i1 = childDesc.hash.i1 + ordValue(childDesc.lastChar) * h;
    childDesc.hash.i2 = childDesc.hash.i2 - (ValueSize<TAlphabet>::VALUE - 1 - ordValue(childDesc.lastChar)) * h;

#ifdef SEQAN_DEBUG
    std::cout << "[" << childDesc.hash.i1 << "," << childDesc.hash.i2 << ")" << std::endl;
#endif

    childDesc.range.i1 = dir[getBucket(index.bucketMap, childDesc.hash.i1)];
    childDesc.range.i2 = dir[getBucket(index.bucketMap, childDesc.hash.i2)];

    childDesc.repLen++;

#ifdef SEQAN_DEBUG
    std::cout << "[" << childDesc.range.i1 << "," << childDesc.range.i2 << ")" << std::endl;
#endif

    return childDesc.range.i1 < childDesc.range.i2;
}

template <typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec>
inline bool _getNextNode(Iter<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >, VSTree<TSpec> > & it)
{
    typedef Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >   TIndex;
    typedef typename Fibre<TIndex, QGramDir>::Type              TDir;
    typedef typename Fibre<TIndex, QGramShape>::Type            TShape;
    typedef typename Host<TShape>::Type                         TAlphabet;
    typedef typename Value<TShape>::Type                        THValue;

    TIndex const & index = container(it);
    TDir const & dir = indexDir(index);
    TShape const & shape = indexShape(index);

    // TODO(esiragusa): Remove call to pow
    THValue h = static_cast<THValue>(pow(static_cast<double>(ValueSize<TAlphabet>::VALUE), static_cast<int>(weight(shape) - value(it).repLen)));

    // TODO(esiragusa): Remove workaround for alphabets with quality values
    for (typename ValueSize<TAlphabet>::Type c = ordValue(value(it).lastChar) + 1; c < ValueSize<TAlphabet>::VALUE; ++c)
    {
        value(it).hash.i1 = value(it).hash.i2;
        value(it).hash.i2 += h;
        value(it).range.i1 = value(it).range.i2;
        value(it).range.i2 = dir[getBucket(index.bucketMap, value(it).hash.i2)];
        value(it).lastChar = c;

        if (value(it).range.i1 < value(it).range.i2)
            return true;
    }

    return false;
}

template <typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec, typename TString, typename TSize>
inline bool
_goDownString(Iter<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > & node,
              TString const & pattern,
              TSize & lcp)
{
    typedef typename Iterator<TString const, Standard>::Type    PatternIterator;

    lcp = 0;

    PatternIterator patternBegin = begin(pattern, Standard());
    PatternIterator patternEnd = end(pattern, Standard());

    // TODO(esiragusa): Specialize to compute the whole hash at once
    for (PatternIterator patternIt = patternBegin; patternIt != patternEnd; ++patternIt)
    {
        if (!_goDownChar(node, *patternIt))
            return false;

        ++lcp;
    }

    return true;
}

template <typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder, typename THideEmptyEdges>
inline bool _goDown(Iter<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > & it,
                    VSTreeIteratorTraits<TDfsOrder, THideEmptyEdges> const)
{
    typedef Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >   TIndex;
    typedef typename Fibre<TIndex, QGramShape>::Type            TShape;
    typedef typename Host<TShape>::Type                         TAlphabet;

    if (isLeaf(it))
        return false;

    // TODO(esiragusa): check nodeHullPredicate

    typename VertexDescriptor<TIndex>::Type nodeDesc;

    // TODO(esiragusa): Remove workaround for alphabets with quality values
    for (typename ValueSize<TAlphabet>::Type c = 0; c < ValueSize<TAlphabet>::VALUE; ++c)
    {
        if (_getNodeByChar(it, c, nodeDesc))
        {
            _historyPush(it);
            value(it) = nodeDesc;
            return true;
        }
    }

    return false;
}

template <typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder, typename THideEmptyEdges>
inline bool _goRight(Iter<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > & it,
                     VSTreeIteratorTraits<TDfsOrder, THideEmptyEdges> const)
{
    if (isRoot(it))
        return false;

    return _getNextNode(it);
}

template <typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec>
inline bool
_goUp(Iter<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > & it)
{
    if (!empty(it.history))
    {
        value(it).range = back(it.history).range;
        value(it).hash = back(it.history).hash;
        value(it).lastChar = back(it.history).lastChar;
        value(it).repLen--;
        eraseBack(it.history);
        return true;
    }
    return false;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INDEX_QGRAM_STREE_H_
