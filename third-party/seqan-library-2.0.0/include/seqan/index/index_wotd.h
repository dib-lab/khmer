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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_INDEX_WOTD_H
#define SEQAN_HEADER_INDEX_WOTD_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
// wotd tree index fibres

/*!
 * @defgroup WOTDIndexFibres WOTD Index Fibres
 * @brief Tag to select a specific fibre (e.g. table, object, ...) of an @link
 *        IndexWotd @endlink index.
 *
 * These tags can be used to get @link Fibre @endlink of an @link IndexWotd @endlink.
 *
 * @see Fibre
 * @see Index#getFibre
 * @see IndexWotd
 *
 * @tag WOTDIndexFibres#WotdDir
 * @brief The child table.
 *
 * @tag WOTDIndexFibres#WotdRawSA
 * @brief The raw suffix array.
 *
 * @tag WOTDIndexFibres#WotdText
 * @brief The original text the index should be based on.
 *
 * @tag WOTDIndexFibres#WotdRawText
 * @brief The raw text the index is really based on.
 *
 * @tag WOTDIndexFibres#WotdSA
 * @brief The suffix array.
 */

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexWotd#indexSA
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>getFibre(.., WotdSA)</tt>.
 *
 * @signature TSa indexSA(index);
 *
 * @param[in] index The @link IndexWotd @endlink object holding the fibre.
 *
 * @return TSa A reference to the @link WOTDIndexFibres#WotdSA @endlink fibre (partially sorted suffix array).
 */

/*!
 * @fn IndexWotd#indexDir
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>getFibre(.., WotdDir())</tt>.
 * @signature TFibre indexDir(index);
 *
 * @param[in] index The @link IndexWotd @endlink object holding the fibre.
 *
 * @return TFibre A reference to the @link WOTDIndexFibres#WotdDir @endlink fibre (tree structure).
 */

/*!
 * @fn IndexWotd#saAt
 * @headerfile <seqan/index.h>
 * @note Advanced functionality, not commonly used.
 * @brief Shortcut for <tt>value(indexSA(..), ..)</tt>.
 *
 * @signature TValue saAt(position, index);
 *
 * @param[in] index The @link IndexWotd @endlink object holding the fibre.
 * @param[in] position A position in the array on which the value should be accessed.
 *
 * @return TValue A reference or proxy to the value in the @link WOTDIndexFibres#WotdSA @endlink fibre.
 *                To be more precise, a reference to a position containing a value of type
 *                @link SAValue @endlink is returned (or a proxy).
 */

/*!
 * @fn IndexWotd#dirAt
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>value(indexDir(index), position)</tt>.
 *
 * @signature TFibre dirAt(position, index);
 *
 * @param[in] index    The @link IndexWotd @endlink object holding the fibre.
 * @param[in] position A position in the array on which the value should be accessed.
 *
 * @return TFibre A reference to the @link WOTDIndexFibres#WotdDir @endlink fibre.
 */

    typedef FibreText        WotdText;
    typedef FibreRawText    WotdRawText;
    typedef FibreSA         WotdSA;
    typedef FibreRawSA        WotdRawSA;
    typedef FibreDir        WotdDir;

//////////////////////////////////////////////////////////////////////////////
// wotd tree index

/*!
 * @class IndexWotd
 * @extends Index
 * @implements StringTreeConcept
 * @headerfile <seqan/index.h>
 * @brief An index based on a lazy suffix tree (see Giegerich et al., "Efficient implementation of lazy suffix
 *        trees").
 *
 * @signature template <typename TText, typename TSpec>
 *            class Index<TText, IndexWotd<TSpec> >;
 *
 * @tparam TText The @link TextConcept @endlink text type.
 * @tparam TSpec The type for further specialization of the Index type.
 *
 * The fibres (see @link Index @endlink and @link Fibre @endlink) of this index are a partially sorted suffix array
 * (see @link WOTDIndexFibres#WotdSA @endlink) and the wotd tree (see @link WOTDIndexFibres#WotdDir @endlink).
 *
 * Demo: Demo.Constraint Iterator
 *
 * @see WOTDIndexFibres
 */

    struct WotdOriginal_;
    typedef Tag<WotdOriginal_> const WotdOriginal;

    template < typename TSpec = void >
    struct IndexWotd {};

/*
    template < typename TObject, typename TSpec >
    struct Fibre< Index<TObject, IndexWotd<TSpec> >, FibreDir>
    {
        typedef Index<TObject, IndexWotd<TSpec> > TIndex;
        typedef String<
            typename typename Size<TIndex>::Type,
            Alloc<>
        > Type;
    };
*/

    template < typename TObject, typename TSpec >
    class Index<TObject, IndexWotd<TSpec> > {
    public:
        typedef typename Fibre<Index, WotdText>::Type        TText;
        typedef typename Fibre<Index, WotdSA>::Type         TSA;
        typedef typename Fibre<Index, WotdDir>::Type        TDir;

        typedef typename Value<Index>::Type                    TValue;
        typedef typename Value<TDir>::Type                    TDirValue;
        typedef typename Size<Index>::Type                    TSize;
        typedef String<TSize, Alloc<> >                        TCounter;
        typedef String<typename Value<TSA>::Type, Alloc<> >    TTempSA;
        typedef typename Cargo<Index>::Type                    TCargo;

        // 1st word flags
        static TDirValue const LEAF          = (TDirValue)1 << (BitsPerValue<TDirValue>::VALUE - 1); // this node is a leaf
        static TDirValue const LAST_CHILD    = (TDirValue)1 << (BitsPerValue<TDirValue>::VALUE - 2); // this node is the last child
        // 2nd word flag
        static TDirValue const UNEVALUATED   = (TDirValue)1 << (BitsPerValue<TDirValue>::VALUE - 1); // this node is partially evalutated and has no evaluated children
        static TDirValue const SENTINELS     = (TDirValue)1 << (BitsPerValue<TDirValue>::VALUE - 2); // the children of this node have solely $-edges

        static TDirValue const BITMASK0      = ~(LEAF | LAST_CHILD);
        static TDirValue const BITMASK1      = ~(UNEVALUATED | SENTINELS);


        Holder<TText>    text;    // underlying text
        TSA                sa;        // suffix array sorted by the first q chars
        TDir            dir;    // bucket directory
        TCargo            cargo;    // user-defined cargo


        TTempSA            tempSA;
        TCounter        tempOcc;
        TCounter        tempBound;

        TSize            sentinelOcc;
        TSize            sentinelBound;
        bool            interSentinelNodes;    // should virtually one (true) $-sign or many (false) $_i-signs be appended to the strings in text

        Index():
            interSentinelNodes(false) {}

        Index(Index &other) :
            text(other.text),
            sa(other.sa),
            dir(other.dir),
            cargo(other.cargo),
            tempSA(other.tempSA),
            tempBound(other.tempBound),
            sentinelOcc(other.sentinelOcc),
            sentinelBound(other.sentinelBound),
            interSentinelNodes(other.interSentinelNodes) {}

        Index(Index const &other) :
            text(other.text),
            sa(other.sa),
            dir(other.dir),
            cargo(other.cargo),
            tempSA(other.tempSA),
            tempBound(other.tempBound),
            sentinelOcc(other.sentinelOcc),
            sentinelBound(other.sentinelBound),
            interSentinelNodes(other.interSentinelNodes) {}

        template <typename TText_>
        Index(TText_ &_text) :
            text(_text),
            sentinelOcc(0),
            sentinelBound(0),
            interSentinelNodes(false) {}

        template <typename TText_>
        Index(TText_ const &_text):
            text(_text),
            sentinelOcc(0),
            sentinelBound(0),
            interSentinelNodes(false) {}
    };
/*
    template < typename TText, typename TSpec >
    struct Value< Index<TText, IndexWotd<TSpec> > > {
        typedef typename Value< typename Fibre< Index<TText, IndexWotd<TSpec> >, WotdRawText >::Type >::Type Type;
    };

    template < typename TText, typename TSpec >
    struct Size< Index<TText, IndexWotd<TSpec> > > {
        typedef typename Size< typename Fibre< Index<TText, IndexWotd<TSpec> >, WotdRawText >::Type >::Type Type;
    };
*/

template <typename TText, typename TSpec>
SEQAN_CONCEPT_IMPL((Index<TText, IndexWotd<TSpec> >), (StringTreeConcept));

template <typename TText, typename TSpec>
SEQAN_CONCEPT_IMPL((Index<TText, IndexWotd<TSpec> > const), (StringTreeConcept));

//////////////////////////////////////////////////////////////////////////////
// default fibre creators

    template < typename TText, typename TSpec >
    struct DefaultIndexCreator<Index<TText, IndexWotd<TSpec> >, FibreDir> {
        typedef Default Type;
    };

    template < typename TText, typename TSSetSpec, typename TSpec >
    struct DefaultIndexCreator<Index<StringSet<TText, TSSetSpec>, IndexWotd<TSpec> >, FibreDir> {
        typedef Default Type;
    };

//////////////////////////////////////////////////////////////////////////////
// default finder

    template < typename TText, typename TSpec >
    struct DefaultFinder< Index<TText, IndexWotd<TSpec> > >
    {
        typedef FinderSTree Type;    // standard wotd finder is tree based search
    };

//////////////////////////////////////////////////////////////////////////////


    template <typename TSize>
    struct VertexWotdOriginal_ {
        TSize        node;            // position of current node entry in directory
        TSize        parentRepLen;    // representative length of parent node
        TSize        edgeLen;        // length of edge above current node

        VertexWotdOriginal_() : node(0), parentRepLen(0), edgeLen(0) {}
        VertexWotdOriginal_(MinimalCtor) : node(0), parentRepLen(0), edgeLen(0)
        {
            _setSizeInval(node);
        }
    };

    template <typename TSize>
    struct VertexWotdModified_ {
        TSize        node;            // position of current node entry in directory
        TSize        parentRepLen;    // representative length of parent node
        TSize        edgeLen;        // length of edge above current node
        Pair<TSize> range;            // current SA interval of hits
        TSize        parentRight;    // right boundary of parent node's range (allows to go right)

        VertexWotdModified_() :
            node(0),
            parentRepLen(0),
            edgeLen(0),
            range(0,0),
            parentRight(0)
        {}
        VertexWotdModified_(MinimalCtor) :
            node(0),
            parentRepLen(0),
            edgeLen(0),
            range(0,0),
            parentRight(0)
        {}
        VertexWotdModified_(Pair<TSize> const &otherRange, TSize otherParentRight) :
            node(0),
            parentRepLen(0),
            edgeLen(0),
            range(otherRange),
            parentRight(otherParentRight)
        {}
    };

//////////////////////////////////////////////////////////////////////////////
    template < typename TText >
    struct VertexDescriptor< Index<TText, IndexWotd<WotdOriginal> > > {
        typedef typename Size< Index<TText, IndexWotd<WotdOriginal> > >::Type TSize;
        typedef VertexWotdOriginal_<TSize> Type;
    };

    template < typename TText, typename TSpec >
    struct VertexDescriptor< Index<TText, IndexWotd<TSpec> > > {
        typedef typename Size< Index<TText, IndexWotd<TSpec> > >::Type TSize;
        typedef VertexWotdModified_<TSize> Type;
    };

    template < typename TText, typename TSpec >
    void _indexRequireTopDownIteration(Index<TText, IndexWotd<TSpec> > &index)
    {
        indexRequire(index, WotdDir());
    }

//////////////////////////////////////////////////////////////////////////////
// history stack functions

    template <typename TSize>
    struct HistoryStackWotdOriginal_
    {
        TSize        node;        // position of current node entry in directory
        TSize        edgeLen;    // length of edge above current node
    };

    template <typename TSize>
    struct HistoryStackWotdModified_
    {
        TSize        node;        // position of current node entry in directory
        TSize        edgeLen;    // length of edge above current node
        Pair<TSize> range;        // current SA interval of hits
    };

    template < typename TText, typename TSpec >
    struct HistoryStackEntry_< Iter< Index<TText, IndexWotd<WotdOriginal> >, VSTree< TopDown< ParentLinks<TSpec> > > > >
    {
        typedef Index<TText, IndexWotd<WotdOriginal> >    TIndex;
        typedef typename Size<TIndex>::Type                TSize;
        typedef HistoryStackWotdOriginal_<TSize>        Type;
    };

    template < typename TText, typename TIndexSpec, typename TSpec >
    struct HistoryStackEntry_< Iter< Index<TText, IndexWotd<TIndexSpec> >, VSTree< TopDown< ParentLinks<TSpec> > > > >
    {
        typedef Index<TText, IndexWotd<TIndexSpec> >    TIndex;
        typedef typename Size<TIndex>::Type                TSize;
        typedef HistoryStackWotdModified_<TSize>        Type;
    };


    template < typename TText, typename TSpec >
    inline void
    _historyPush(Iter< Index<TText, IndexWotd<WotdOriginal> >, VSTree< TopDown<TSpec> > > &it)
    {
        it._parentDesc = value(it);
        value(it).parentRepLen += parentEdgeLength(it);
    }

    template < typename TText, typename TIndexSpec, typename TSpec >
    inline void
    _historyPush(Iter< Index<TText, IndexWotd<TIndexSpec> >, VSTree< TopDown<TSpec> > > &it)
    {
        it._parentDesc = value(it);
        value(it).parentRepLen += parentEdgeLength(it);
        value(it).parentRight = value(it).range.i2;
    }

    template < typename TText, typename TSpec >
    inline void
    _historyPush(Iter< Index<TText, IndexWotd<WotdOriginal> >, VSTree< TopDown< ParentLinks<TSpec> > > > &it)
    {
        typedef typename Size< Index<TText, IndexWotd<WotdOriginal> > >::Type TSize;
        TSize edgeLen = parentEdgeLength(it);
        HistoryStackWotdOriginal_<TSize> entry = { value(it).node, edgeLen };
        appendValue(it.history, entry);
        value(it).parentRepLen += edgeLen;
    }

    template < typename TText, typename TIndexSpec, typename TSpec >
    inline void
    _historyPush(Iter< Index<TText, IndexWotd<TIndexSpec> >, VSTree< TopDown< ParentLinks<TSpec> > > > &it)
    {
        typedef typename Size< Index<TText, IndexWotd<TIndexSpec> > >::Type TSize;
        TSize edgeLen = parentEdgeLength(it);
        HistoryStackWotdModified_<TSize> entry = { value(it).node, edgeLen, value(it).range };
        appendValue(it.history, entry);
        value(it).parentRepLen += edgeLen;
        value(it).parentRight = value(it).range.i2;
    }

//////////////////////////////////////////////////////////////////////////////
    template < typename TText, typename TIndexSpec, typename TPropertyMap >
    inline void
    resizeVertexMap(
        TPropertyMap & pm,
        Index<TText, IndexWotd<TIndexSpec> > const& index)
    {
        resize(pm, length(indexDir(index)), Generous());
    }

/* // different interface compared to resizeVertexMap(graph, ...)
    template < typename TText, typename TIndexSpec, typename TPropertyMap, typename TProperty >
    inline void
    resizeVertexMap(
        Index<TText, IndexWotd<TIndexSpec> > const& index,
        TPropertyMap & pm,
        TProperty const & prop)
    {
        resize(pm, length(indexDir(index)), prop, Generous());
    }
*/
    template < typename TSize >
    inline typename Id< VertexWotdOriginal_<TSize> const >::Type
    _getId(VertexWotdOriginal_<TSize> const &desc)
    {
        return desc.node;
    }

    template < typename TSize >
    inline typename Id< VertexWotdOriginal_<TSize> >::Type
    _getId(VertexWotdOriginal_<TSize> &desc)
    {
        return _getId(const_cast<VertexWotdOriginal_<TSize> const &>(desc));
    }

    template < typename TSize >
    inline typename Id< VertexWotdModified_<TSize> const >::Type
    _getId(VertexWotdModified_<TSize> const &desc)
    {
        return desc.node;
    }

    template < typename TSize >
    inline typename Id< VertexWotdModified_<TSize> >::Type
    _getId(VertexWotdModified_<TSize> &desc)
    {
        return _getId(const_cast<VertexWotdModified_<TSize> const &>(desc));
    }

//////////////////////////////////////////////////////////////////////////////

    template < typename TSize >
    inline bool _isRoot(VertexWotdOriginal_<TSize> const &value)
    {
        return value.node == 0;
    }

    template < typename TSize >
    inline bool _isRoot(VertexWotdModified_<TSize> const &value)
    {
        return value.node == 0;
    }

    // is this a leaf? (including empty $-edges)
    template < typename TText, typename TIndexSpec, typename TSpec, typename TDfsOrder >
    inline bool _isLeaf(
        Iter< Index<TText, IndexWotd<TIndexSpec> >, VSTree<TSpec> > const &it,
        VSTreeIteratorTraits<TDfsOrder, False> const)
    {
        typedef Index<TText, IndexWotd<TIndexSpec> > TIndex;
        TIndex const &index = container(it);
        return (dirAt(value(it).node, index) & index.LEAF) != 0;
    }

    // is this a leaf? (excluding empty $-edges)
    template < typename TText, typename TIndexSpec, typename TSpec, typename TDfsOrder >
    inline bool _isLeaf(
        Iter< Index<TText, IndexWotd<TIndexSpec> >, VSTree<TSpec> > const &it,
        VSTreeIteratorTraits<TDfsOrder, True> const)
    {
        typedef Index<TText, IndexWotd<TIndexSpec> >    TIndex;

        TIndex const &index = container(it);
        if (dirAt(value(it).node, index) & index.LEAF)
            return true;

        // ensure node evaluation and test for sentinel child edges?
        return _wotdEvaluate(it) & index.SENTINELS;
    }


    // parentEdgeLength - ORIGINAL VERSION
    template < typename TIndex, typename TSize >
    inline typename Size<TIndex>::Type
    parentEdgeLength(TIndex const &index, VertexWotdOriginal_<TSize> &vDesc)
    {
        TSize edgeLen = vDesc.edgeLen;
        if (edgeLen != (TSize)-1)
            return edgeLen;

        TSize pos = vDesc.node;
        TSize w0 = dirAt(pos, index);
        if (w0 & index.LEAF)
            return vDesc.edgeLen = suffixLength(w0 & index.BITMASK0, index);

        TSize w1 = dirAt(pos + 1, index);
        if (w1 & index.UNEVALUATED)
            return vDesc.edgeLen = _bucketLcp(
                infix(indexSA(index), w0 & index.BITMASK0, w1 & index.BITMASK1),
                indexText(index));
        else
            return vDesc.edgeLen = _getNodeLP(index, w1) - (w0 & index.BITMASK0);
    }

    // parentEdgeLength - MODIFIED VERSION
    template < typename TIndex, typename TSize >
    inline typename Size<TIndex>::Type
    parentEdgeLength(TIndex const &index, VertexWotdModified_<TSize> &vDesc)
    {
        TSize edgeLen = vDesc.edgeLen;
        if (edgeLen != (TSize)-1)
            return edgeLen;

        TSize pos = vDesc.node;
        TSize w0 = dirAt(pos, index);
        if (w0 & index.LEAF)
            return vDesc.edgeLen =
                suffixLength(saAt(vDesc.range.i1, index), index) - vDesc.parentRepLen;

        TSize w1 = dirAt(pos + 1, index);
        if (w1 & index.UNEVALUATED)
            if (_isSizeInval(vDesc.range.i2))
                return vDesc.edgeLen = _bucketLcp(
                    suffix(indexSA(index), vDesc.range.i1),
                    indexText(index),
                    vDesc.parentRepLen) - vDesc.parentRepLen;
            else
                return vDesc.edgeLen = _bucketLcp(
                    infix(indexSA(index), vDesc.range.i1, vDesc.range.i2),
                    indexText(index),
                    vDesc.parentRepLen) - vDesc.parentRepLen;
        else
            return (dirAt(w1 & index.BITMASK1, index) & index.BITMASK0) - vDesc.parentRepLen;
    }

    template < typename TText, typename TIndexSpec, typename TSpec >
    inline typename Size< Index<TText, IndexWotd<TIndexSpec> > >::Type
    parentEdgeLength(Iter<
        Index<TText, IndexWotd<TIndexSpec> >,
        VSTree< TopDown<TSpec> > > const &it)
    {
        typedef Iter< Index<TText, IndexWotd<TIndexSpec> >, VSTree< TopDown<TSpec> > > TIter;
        return parentEdgeLength(container(it), value(const_cast<TIter&>(it)));
    }

    template < typename TText, typename TIndexSpec, typename TSpec >
    inline typename Size< Index<TText, IndexWotd<TIndexSpec> > >::Type
    parentRepLength(Iter<
        Index<TText, IndexWotd<TIndexSpec> >,
        VSTree< TopDown<TSpec> > > const &it)
    {
        return value(it).parentRepLen;
    }

    template < typename TText, typename TIndexSpec, typename TSpec >
    inline typename Size< Index<TText, IndexWotd<TIndexSpec> > >::Type
    parentRepLength(Iter<
        Index<TText, IndexWotd<TIndexSpec> >,
        VSTree< TopDown< ParentLinks<TSpec> > > > const &it)
    {
        return value(it).parentRepLen;
    }

    template < typename TText, typename TIndexSpec, typename TSpec >
    inline typename Size< Index<TText, IndexWotd<TIndexSpec> > >::Type
    repLength(Iter<
        Index<TText, IndexWotd<TIndexSpec> >,
        VSTree< TopDown<TSpec> > > const &it)
    {
        return parentRepLength(it) + parentEdgeLength(it);
    }


    // parentEdgeLabel - ORIGINAL VERSION (doesn't work if TText is a StringSet)
    template < typename TText, typename TSpec >
    inline typename Infix< typename Fibre<Index<TText, IndexWotd<WotdOriginal> >, EsaText>::Type const >::Type
    parentEdgeLabel(Iter< Index<TText, IndexWotd<WotdOriginal> >, VSTree< TopDown<TSpec> > > const &it)
    {
        typedef Index<TText, IndexWotd<WotdOriginal> >    TIndex;
        typedef typename Size<TIndex>::Type                TSize;

        TIndex const &index = container(it);

        if (isRoot(it))
            return infix(indexText(index), 0, 0);
        else {
            TSize occ = _getNodeLP(index, value(it).node);
            return infix(indexText(index), occ, occ + parentEdgeLength(it));
        }
    }

    // getOccurrence - ORIGINAL VERSION
    template < typename TText, typename TSpec >
    inline typename SAValue<Index<TText, IndexWotd<WotdOriginal> > >::Type
    getOccurrence(Iter< Index<TText, IndexWotd<WotdOriginal> >, VSTree<TSpec> > const &it) {
        return _getNodeLP(container(it), value(it).node) - value(it).parentRepLen;
    }

    template < typename TText, typename TIndexSpec, typename TSpec >
    inline bool
    emptyParentEdge(Iter< Index<TText, IndexWotd<TIndexSpec> >, VSTree<TopDown<TSpec> > > const &it)
    {
        typedef Index<TText, IndexWotd<TIndexSpec> > TIndex;

        TIndex const &index = container(it);
        typename SAValue<TIndex>::Type pos = getOccurrence(it);
        return getSeqOffset(pos, stringSetLimits(index)) + value(it).parentRepLen
            == sequenceLength(getSeqNo(pos, stringSetLimits(index)), index);
    }

    // to avoid ambiguity
    template < typename TText, typename TIndexSpec, typename TSpec >
    inline bool
    emptyParentEdge(Iter< Index<TText, IndexWotd<TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > const &it)
    {
        typedef Index<TText, IndexWotd<TIndexSpec> > TIndex;

        TIndex const &index = container(it);
        typename SAValue<TIndex>::Type pos = getOccurrence(it);
        return getSeqOffset(pos, stringSetLimits(index)) + value(it).parentRepLen
            == sequenceLength(getSeqNo(pos, stringSetLimits(index)), index);
    }



    template < typename TText, typename TSpec >
    inline void
    goRoot(Iter<
        Index<TText, IndexWotd<WotdOriginal> >,
        VSTree<TSpec> > &it)
    {
        _historyClear(it);
        value(it).node = 0;            // start in root node (first entry in dir)
        value(it).parentRepLen = 0;    // parent prefix length is 0
        value(it).edgeLen = 0;        // edge length is 0
    }

    template < typename TText, typename TIndexSpec, typename TSpec >
    inline void
    goRoot(Iter<
        Index<TText, IndexWotd<TIndexSpec> >,
        VSTree<TSpec> > &it)
    {
        _historyClear(it);
        value(it).range.i1 = 0;        // start in root node with range (0,infty)
        _setSizeInval(value(it).range.i2);    // infty is equivalent to length(index) and faster to compare
        value(it).node = 0;            // start in root node (first entry in dir)
        value(it).parentRepLen = 0;    // parent prefix length is 0
        value(it).edgeLen = 0;        // edge length is 0
    }

    template < typename TText, typename TSpec >
    inline bool atEnd(Iter<Index<TText, IndexWotd<WotdOriginal> >, VSTree<TSpec> > &it)
    {
        return _isSizeInval(value(it).node);
    }

    template < typename TText, typename TSpec >
    inline bool atEnd(Iter<Index<TText, IndexWotd<WotdOriginal> >, VSTree<TSpec> > const &it)
    {
        return _isSizeInval(value(it).node);
    }


    // adjust iterator's right border of SA range
    template < typename TText, typename TSpec >
    inline void
    _adjustRightBorder(
        Iter< Index<TText, IndexWotd<WotdOriginal> >, VSTree< TopDown<TSpec> > > &)
    {}

    template < typename TText, typename TIndexSpec, typename TSpec >
    inline void
    _adjustRightBorder(
        Iter< Index<TText, IndexWotd<TIndexSpec> >, VSTree< TopDown<TSpec> > > &it)
    {
        typedef Index<TText, IndexWotd<TIndexSpec> >    TIndex;
        typedef typename Size<TIndex>::Type                TSize;

        TIndex    const &index = container(it);
        TSize    pos = value(it).node;

        TSize w0 = dirAt(pos, index);
        if (w0 & index.LEAF)
            value(it).range.i2 = value(it).range.i1 + 1;
        else
            if (w0 & index.LAST_CHILD)
                value(it).range.i2 = value(it).parentRight;
            else {
                w0 = dirAt(pos + 2, index);
                value(it).range.i2 = w0 & index.BITMASK0;
            }
    }

    // go down the leftmost edge (including empty $-edges)
    template < typename TText, typename TSpec, typename TDfsOrder, typename THideEmptyEdges >
    inline bool
    _goDown(
        Iter< Index<TText, IndexWotd<WotdOriginal> >, VSTree< TopDown<TSpec> > > &it,
        VSTreeIteratorTraits<TDfsOrder, THideEmptyEdges> const)
    {
        typedef Index<TText, IndexWotd<WotdOriginal> >    TIndex;
        typedef typename Size<TIndex>::Type                TSize;

        if (_isLeaf(it, EmptyEdges())) return false;

        TIndex &index = container(it);
        _historyPush(it);

        // ensure node evaluation
        TSize childNode = _wotdEvaluate(it);

        if (THideEmptyEdges::VALUE && (childNode & index.SENTINELS) != 0)
            return false;

        // go down
        value(it).node = childNode & index.BITMASK1;
        value(it).edgeLen = -1;

        // go right if parent edge is empty
        // or hull predicate is false
        if ((THideEmptyEdges::VALUE && emptyParentEdge(it)) || !nodeHullPredicate(it))
            if (!goRight(it)) {
                _goUp(it);
                return false;
            }

        return true;
    }

    // go down the leftmost edge (excluding empty $-edges)
    template < typename TText, typename TIndexSpec, typename TSpec, typename TDfsOrder, typename THideEmptyEdges >
    inline bool
    _goDown(
        Iter< Index<TText, IndexWotd<TIndexSpec> >, VSTree< TopDown<TSpec> > > &it,
        VSTreeIteratorTraits<TDfsOrder, THideEmptyEdges> const)
    {
        typedef Index<TText, IndexWotd<TIndexSpec> >    TIndex;
        typedef typename Size<TIndex>::Type                TSize;

        if (_isLeaf(it, EmptyEdges())) return false;
        TIndex const &index = container(it);

        // ensure node evaluation
        TSize childNode = _wotdEvaluate(it);

        if (THideEmptyEdges::VALUE && (childNode & index.SENTINELS) != 0)
            return false;

        // go down
        _historyPush(it);
        value(it).node = childNode & index.BITMASK1;
        value(it).edgeLen = -1;
        _adjustRightBorder(it);

        // go right if parent edge is empty
        // or hull predicate is false
        if ((THideEmptyEdges::VALUE && emptyParentEdge(it)) || !nodeHullPredicate(it))
            if (!goRight(it)) {
                _goUp(it);
                return false;
            }

        return true;
    }

    // go right to the lexic. next sibling
    template < typename TText, typename TSpec, typename TDfsOrder, typename THideEmptyEdges >
    inline bool
    _goRight(
        Iter< Index<TText, IndexWotd<WotdOriginal> >, VSTree< TopDown<TSpec> > > &it,
        VSTreeIteratorTraits<TDfsOrder, THideEmptyEdges> const)
    {
        typedef Index<TText, IndexWotd<WotdOriginal> >    TIndex;
        typedef typename Size<TIndex>::Type                TSize;

        TIndex const &index = container(it);

        do {
            TSize w0 = dirAt(value(it).node, index);
            if (w0 & index.LAST_CHILD)
                return false;

            value(it).node += (w0 & index.LEAF)? 1: 2;
            value(it).edgeLen = -1;

            _adjustRightBorder(it);
        } while ((THideEmptyEdges::VALUE && emptyParentEdge(it)) || !nodeHullPredicate(it));
        return true;
    }

    template < typename TText, typename TIndexSpec, typename TSpec, typename TDfsOrder, typename THideEmptyEdges >
    inline bool
    _goRight(
        Iter< Index<TText, IndexWotd<TIndexSpec> >, VSTree< TopDown<TSpec> > > &it,
        VSTreeIteratorTraits<TDfsOrder, THideEmptyEdges> const)
    {
        typedef Index<TText, IndexWotd<TIndexSpec> >    TIndex;
        typedef typename Size<TIndex>::Type                TSize;

        TIndex const &index = container(it);

        do {
            TSize w0 = dirAt(value(it).node, index);
            if (w0 & index.LAST_CHILD)
                return false;

            value(it).node += (w0 & index.LEAF)? 1: 2;
            value(it).edgeLen = -1;

            value(it).range.i1 = value(it).range.i2;
            _adjustRightBorder(it);
        } while ((THideEmptyEdges::VALUE && emptyParentEdge(it)) || !nodeHullPredicate(it));
        return true;
    }

    // go up one edge (returns false if in root node)
    // can be used at most once, as no history stack is available
    template < typename TText, typename TWotdSpec, typename TSpec >
    inline bool
    _goUp(Iter< Index<TText, IndexWotd<TWotdSpec> >, VSTree< TopDown<TSpec> > > &it)
    {
        if (!isRoot(it)) {
            value(it) = it._parentDesc;
            return true;
        }
        return false;
    }

    // go up one edge (returns false if in root node)
    template < typename TText, typename TSpec >
    inline bool
    _goUp(Iter< Index<TText, IndexWotd<WotdOriginal> >, VSTree< TopDown< ParentLinks<TSpec> > > > &it)
    {
        typedef typename Size< Index<TText, IndexWotd<WotdOriginal> > >::Type TSize;

        if (!empty(it.history)) {
            HistoryStackWotdOriginal_<TSize> const &entry = back(it.history);
            value(it).node = entry.node;
            value(it).parentRepLen -= entry.edgeLen;
            value(it).edgeLen = entry.edgeLen;
            eraseBack(it.history);
            return true;
        }
        return false;
    }

    template < typename TText, typename TIndexSpec, typename TSpec >
    inline bool
    _goUp(Iter< Index<TText, IndexWotd<TIndexSpec> >, VSTree< TopDown< ParentLinks<TSpec> > > > &it)
    {
        typedef typename Size< Index<TText, IndexWotd<TIndexSpec> > >::Type TSize;

        if (!empty(it.history))
        {
            HistoryStackWotdModified_<TSize> const &entry = back(it.history);
            value(it).node = entry.node;
            value(it).parentRepLen -= entry.edgeLen;
            value(it).edgeLen = entry.edgeLen;
            value(it).range = entry.range;
            eraseBack(it.history);
            if (!empty(it.history))
                value(it).parentRight = back(it.history).range.i2;    // copy right boundary of parent's range
            return true;
        }
        return false;
    }

    // return vertex descriptor of parent's node
    template < typename TText, typename TSpec >
    inline typename VertexDescriptor< Index<TText, IndexWotd<WotdOriginal> > >::Type
    nodeUp(Iter< Index<TText, IndexWotd<WotdOriginal> >, VSTree< TopDown< ParentLinks<TSpec> > > > const &it)
    {
        typedef Index<TText, IndexWotd<WotdOriginal> > TIndex;
        typedef typename Size<TIndex>::Type TSize;

        if (!empty(it.history))
        {
            HistoryStackWotdOriginal_<TSize> const &entry = back(it.history);
            typename VertexDescriptor<TIndex>::Type desc;

            desc.node = entry.node;
            desc.parentRepLen = value(it).parentRepLen - entry.edgeLen;
            desc.edgeLen = entry.edgeLen;
            TSize h = length(it.history) - 1;
            if (h != 0) --h;
            return desc;
        } else
            return value(it);
    }

    // return vertex descriptor of parent's node
    template < typename TText, typename TIndexSpec, typename TSpec >
    inline typename VertexDescriptor< Index<TText, IndexWotd<TIndexSpec> > >::Type
    nodeUp(Iter< Index<TText, IndexWotd<TIndexSpec> >, VSTree< TopDown< ParentLinks<TSpec> > > > const &it)
    {
        typedef Index<TText, IndexWotd<TIndexSpec> > TIndex;
        typedef typename Size<TIndex>::Type TSize;

        if (!empty(it.history))
        {
            HistoryStackWotdModified_<TSize> const &entry = back(it.history);
            typename VertexDescriptor<TIndex>::Type desc;

            desc.node = entry.node;
            desc.parentRepLen = value(it).parentRepLen - entry.edgeLen;
            desc.edgeLen = entry.edgeLen;
            return desc;
        } else
            return value(it);
    }

//////////////////////////////////////////////////////////////////////////////


    // Counting sort - Step 2a: Count characters
    template < typename TBuckets, typename TText >
    inline void
    _wotdCountChars(TBuckets &buckets, TText const &text)
    {
        typedef typename Iterator<TText const, Standard>::Type TTextIterator;

        TTextIterator itText = begin(text, Standard());
        TTextIterator itTextEnd = end(text, Standard());
        for (; itText != itTextEnd; ++itText)
            ++buckets[ordValue(getValue(itText))];
    }


    template < typename TBuckets, typename TText, typename TSpec >
    inline void
    _wotdCountChars(TBuckets &buckets, StringSet<TText, TSpec> const &stringSet)
    {
        typedef typename Iterator<TText const, Standard>::Type TTextIterator;

        for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo)
        {
            TText const &text = value(stringSet, seqNo);
            TTextIterator itText = begin(text, Standard());
            TTextIterator itTextEnd = end(text, Standard());
            for (; itText != itTextEnd; ++itText)
                ++buckets[ordValue(getValue(itText))];
        }
    }

    // Counting sort - Step 2b: Count the prefixLen'th characters of suffices
    template < typename TBuckets, typename TText, typename TSA, typename TSize >
    inline typename Size<TText>::Type
    _wotdCountCharsWotdOriginal(
        TBuckets &buckets,
        TText const &text,
        TSA &sa,
        TSize prefixLen)
    {
        typedef typename Iterator<TText const, Standard>::Type    TTextIterator;
        typedef typename Iterator<TSA, Standard>::Type          TSAIterator;
        typedef typename Size<TText>::Type                        TTextSize;

        TTextIterator itText = begin(text, Standard());
        TSAIterator itSA = begin(sa, Standard());
        TSAIterator itSAEnd = end(sa, Standard());

        TTextSize sentinels = 0;
        TTextSize textLength = length(text);
        for (; itSA != itSAEnd; ++itSA)
        {
            // add prefix on entries in sa
            TTextSize saValue = (*itSA += prefixLen);
            if (textLength > saValue)
                ++buckets[ordValue(*(itText + saValue))];
            else
                if (textLength == saValue) ++sentinels;
        }
        return sentinels;
    }

    template < typename TBuckets, typename TText, typename TSA, typename TSize >
    inline typename Size<TText>::Type
    _wotdCountChars(
        TBuckets &buckets,
        TText const &text,
        TSA const &sa,
        TSize prefixLen)
    {
        typedef typename Iterator<TText const, Standard>::Type    TTextIterator;
        typedef typename Iterator<TSA const, Standard>::Type    TSAIterator;
        typedef typename Size<TText>::Type                        TTextSize;

        TTextIterator itText = begin(text, Standard()) + prefixLen;
        TSAIterator itSA = begin(sa, Standard());
        TSAIterator itSAEnd = end(sa, Standard());

        TTextSize sentinels = 0;
        TTextSize textLength = length(text) - prefixLen;
        for (; itSA != itSAEnd; ++itSA)
        {
            TTextSize saValue = *itSA;
            if (textLength > saValue)
                ++buckets[ordValue(*(itText + saValue))];
            else
                if (textLength == saValue) ++sentinels;
        }
        return sentinels;
    }

    template <
        typename TBuckets,
        typename TText,
        typename TSpec,
        typename TSA,
        typename TSize
    >
    inline typename Size<TText>::Type
    _wotdCountChars(
        TBuckets &buckets,
        StringSet<TText, TSpec> const &stringSet,
        TSA const &sa,
        TSize prefixLen)
    {
        typedef typename Iterator<TText const, Standard>::Type    TTextIterator;
        typedef typename Iterator<TSA const, Standard>::Type    TSAIterator;
        typedef typename Size<TText>::Type                        TTextSize;

        if (empty(stringSet))
            return 0;

        TTextIterator itText = begin(front(stringSet), Standard());
        TSAIterator itSA = begin(sa, Standard());
        TSAIterator itSAEnd = end(sa, Standard());

        TTextSize sentinels = 0;
        TTextSize textLength = 0;
        unsigned lastSeqSeen = (unsigned)-1;
        Pair<unsigned, TTextSize> lPos;
        for (; itSA != itSAEnd; ++itSA)
        {
            posLocalize(lPos, *itSA, stringSetLimits(stringSet));
            if (lastSeqSeen != getSeqNo(lPos))
            {
                lastSeqSeen = getSeqNo(lPos);

                // shift textBegin and textLength by prefixLen
                textLength = length(stringSet[lastSeqSeen]) - prefixLen;
                itText = begin(stringSet[lastSeqSeen], Standard()) + prefixLen;
            }
            if (textLength > getSeqOffset(lPos))
                ++buckets[ordValue(*(itText + getSeqOffset(lPos)))];
            else
                if (textLength == getSeqOffset(lPos)) ++sentinels;
        }
        return sentinels;
    }


//////////////////////////////////////////////////////////////////////////////


    // Counting sort - Step 3: Cumulative sum
    template < typename TBounds, typename TBuckets, typename TSize >
    inline typename Size<TBuckets>::Type
    _wotdCummulativeSum(TBounds &bounds, TBuckets const &buckets, TSize offset)
    {
        typedef typename Iterator<TBounds, Standard>::Type          TBoundIterator;
        typedef typename Iterator<TBuckets const, Standard>::Type   TBucketsIterator;

        TBucketsIterator it = begin(buckets, Standard());
        TBucketsIterator itEnd = end(buckets, Standard());
        TBoundIterator bit = begin(bounds, Standard());

        typename Value<TBounds>::Type    sum = offset;
        typename Size<TBuckets>::Type    requiredSize = 0;
        typename Value<TBuckets>::Type    diff;

        for (; it != itEnd; ++it, ++bit)
            if ((diff = *it) != 0) {
                requiredSize += (diff > 1)? 2: 1;
                *bit = sum;
                sum += diff;
            }

        return requiredSize;
    }

    template < typename TBounds, typename TBuckets >
    inline typename Size<TBuckets>::Type
    _wotdCummulativeSum(TBounds &bounds, TBuckets const &buckets)
    {
        return _wotdCummulativeSum(bounds, buckets, 0);
    }

//////////////////////////////////////////////////////////////////////////////
//TODO(singer): The function createWotdIndex in never defined!
/*!
 * @fn IndexWotd#createWotdIndex
 * @headerfile <seqan/index.h>
 * @brief Builds a the WOTD index.
 *
 * @signature void createWotdIndex(sa, dir, text);
 *
 * @param[out] sa  The resulting list in which all <i>q</i>-grams are sorted alphabetically.
 * @param[out] dir The resulting array that indicates at which position in index the corresponding <i>q</i>-grams
 *                 can be found.
 * @param[in] text The sequence. Types: @link ContainerConcept @endlink
 *
 * The resulting <tt>index</tt> contains the sorted list of qgrams.  For each possible <i>q</i>-gram pos contains
 * the first position in index that corresponds to this <i>q</i>-gram.
 */

    // single sequence
    template < typename TIndex >
    typename Size<TIndex>::Type
    _sortFirstWotdBucket(TIndex &index)
    {
        typedef typename Fibre<TIndex, WotdText >::Type        TText;
        typedef typename Fibre<TIndex, WotdSA >::Type            TSA;
        typedef typename TIndex::TCounter                        TCounter;

        typedef typename Iterator<TText const, Standard>::Type    TTextIterator;
        typedef typename Iterator<TSA, Standard>::Type            TSAIterator;
        typedef typename Iterator<TCounter, Standard>::Type        TCntIterator;
        typedef typename Size<TText>::Type                        TSize;

        TText const &text = indexText(index);
        TCounter &occ = index.tempOcc;
        TCounter &bound = index.tempBound;

        // 1. clear counters and copy SA to temporary SA
        arrayFill(begin(occ, Standard()), end(occ, Standard()), 0);

        // 2. count characters
        _wotdCountChars(occ, text);

        // 3. cumulative sum
        TSize requiredSize = _wotdCummulativeSum(bound, occ);

        // 4. fill suffix array
        {
            TSA &sa = indexSA(index);
            TSAIterator saBeg = begin(sa, Standard());
            TCntIterator boundBeg = begin(bound, Standard());

            TTextIterator itText = begin(text, Standard());
            TTextIterator itTextEnd = end(text, Standard());
            for(TSize i = 0; itText != itTextEnd; ++itText, ++i)
                *(saBeg + (*(boundBeg + ordValue(getValue(itText))))++) = i;
        }
        index.sentinelOcc = 0;
        index.sentinelBound = 0;

        return requiredSize;
    }

    // multiple sequences
    template < typename TText, typename TSpec, typename TIndexSpec >
    typename Size< Index<StringSet<TText, TSpec>, TIndexSpec> >::Type
    _sortFirstWotdBucket(Index<StringSet<TText, TSpec>, TIndexSpec> &index)
    {
        typedef Index<StringSet<TText, TSpec>, TIndexSpec>        TIndex;
        typedef typename Fibre<TIndex, WotdSA >::Type            TSA;
        typedef typename TIndex::TCounter                        TCounter;

        typedef typename Iterator<TText const, Standard>::Type    TTextIterator;
        typedef typename Iterator<TSA, Standard>::Type            TSAIterator;
        typedef typename Iterator<TCounter, Standard>::Type        TCntIterator;
        typedef typename Size<TText>::Type                        TSize;

        StringSet<TText, TSpec> const &stringSet = indexText(index);
        TCounter &occ = index.tempOcc;
        TCounter &bound = index.tempBound;

        // 1. clear counters and copy SA to temporary SA
        arrayFill(begin(occ, Standard()), end(occ, Standard()), 0);

        // 2. count characters
        _wotdCountChars(occ, stringSet);

        // 3. cummulative sum
        TSize requiredSize = _wotdCummulativeSum(bound, occ);

        // 4. fill suffix array
        for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo)
        {
            TSA &sa = indexSA(index);
            TSAIterator saBeg = begin(sa, Standard());
            TCntIterator boundBeg = begin(bound, Standard());

            typename Value<TSA>::Type localPos;
            assignValueI1(localPos, seqNo);
            assignValueI2(localPos, 0);

            TText const &text = value(stringSet, seqNo);
            TTextIterator itText = begin(text, Standard());
            TTextIterator itTextEnd = end(text, Standard());
            for(; itText != itTextEnd; ++itText)
            {
                *(saBeg + (*(boundBeg + ordValue(getValue(itText))))++) = localPos;
                assignValueI2(localPos, getValueI2(localPos) + 1);
            }
        }
        index.sentinelOcc = 0;
        index.sentinelBound = 0;

        return requiredSize;
    }



    // sort bucket using counting sort
    // (nearly) ORIGINAL VERSION:
    // - brings the bucket with the longest suffix (lpBucket) to front
    // - all other buckets are in lexicographical order
    // - adds prefixLen to all SA entries in SA[left,right)
    template < typename TText, typename TBeginPos, typename TEndPos, typename TSize >
    TSize _sortWotdBucket(
        Index<TText, IndexWotd<WotdOriginal> > &index,
        TBeginPos left,
        TEndPos right,
        TSize prefixLen)
    {
        typedef Index<TText, IndexWotd<WotdOriginal> >              TIndex;
        typedef typename Fibre<TIndex, WotdSA >::Type               TSA;
        typedef typename TIndex::TCounter                           TCounter;

        typedef typename Iterator<TText const, Standard>::Type      TTextIterator;
        typedef typename Iterator<TSA, Standard>::Type              TSAIterator;
        typedef typename Iterator<TCounter, Standard>::Type         TCntIterator;
        typedef typename Iterator<TCounter const, Standard>::Type   TConstCntIterator;
        typedef typename Size<TText>::Type                          TTextSize;

        TText const &text = indexText(index);
        TCounter const &tempSA = index.tempSA;
        TCounter &occ = index.tempOcc;
        TCounter &bound = index.tempBound;

        // 1. clear counters and copy SA to temporary SA
        arrayFill(begin(occ, Standard()), end(occ, Standard()), 0);

        // 2. count characters
        index.tempSA = infix(indexSA(index), left, right);
        index.sentinelBound = 0;
        index.sentinelOcc =
            _wotdCountCharsWotdOriginal(occ, text, index.tempSA, prefixLen);

        // 3. cumulative sum
        TSize requiredSize = 0;

        // actually, here sentinelOcc<=1 holds (this is the original wotd)
        if (index.interSentinelNodes) {
            if (index.sentinelOcc != 0)
                requiredSize = (index.sentinelOcc > 1)? 2: 1;    // insert *one* $-edge for all $_i suffices
        } else
            requiredSize = index.sentinelOcc;                    // insert each $_i suffix one-by-one

        requiredSize += _wotdCummulativeSum(bound, occ, left + index.sentinelOcc);
        index.sentinelBound = left;

        // 4. fill suffix array
        {
            TSA &sa = indexSA(index);
            TSAIterator saBeg = begin(sa, Standard());
            TCntIterator boundBeg = begin(bound, Standard());

            TTextIterator itText = begin(text, Standard());
            TConstCntIterator itSA = begin(tempSA, Standard());
            TConstCntIterator itSAEnd = end(tempSA, Standard());
            TTextSize textLength = length(text);
            for(; itSA != itSAEnd; ++itSA)
            {
                TTextSize saValue = *itSA;
                if (textLength > saValue)
                    *(saBeg + (*(boundBeg + ordValue(*(itText + saValue))))++) = saValue;
                else
                    if (textLength == saValue)
                        *(saBeg + index.sentinelBound++) = saValue;
            }
        }

        // 5. move lpBucket to front and add saOffset to hash directory entries
        {
            TSize lpBucket = ordValue(text[tempSA[0]]);
            if (lpBucket != 0) {
                TSize lpBucketOcc = occ[lpBucket];
                TSize lpBucketBound = bound[lpBucket];

                TCntIterator itOcc = begin(occ, Standard()) + lpBucket;
                TCntIterator itBound = begin(bound, Standard()) + lpBucket;
                TCntIterator itBeg = begin(bound, Standard());
                for(; itBound != itBeg; --itBound, --itOcc) {
                    *itOcc = *(itOcc - 1);
                    *itBound = *(itBound - 1);
                }
                if (index.sentinelOcc != 0) {
                    // bring first bucket before sentinel bucket
                    *itOcc = index.sentinelOcc;
                    *itBound = index.sentinelBound;
                    index.sentinelOcc = lpBucketOcc;
                    index.sentinelBound = lpBucketBound;
                } else {
                    *itOcc = lpBucketOcc;
                    *itBound = lpBucketBound;
                }
            } else
                if (index.sentinelOcc != 0) {
                    // bring first bucket before sentinel bucket
                    TSize swap = index.sentinelOcc;
                    index.sentinelOcc = occ[0];
                    occ[0] = swap;
                    swap = index.sentinelBound;
                    index.sentinelBound = bound[0];
                    bound[0] = swap;
                }
        }

        return requiredSize;
    }





    // sort bucket using counting sort
    // MODIFIED VERSION:
    // - all buckets are in lexicographical order
    // - SA[left,right) contains real SA entries (the beginning positions of the suffices)

    // single sequence
    template < typename TIndex, typename TBeginPos, typename TEndPos, typename TSize >
    TSize _sortWotdBucket(
        TIndex &index,
        TBeginPos left,
        TEndPos right,
        TSize prefixLen)
    {
        typedef typename Fibre<TIndex, WotdText >::Type             TText;
        typedef typename Fibre<TIndex, WotdSA >::Type               TSA;
        typedef typename TIndex::TCounter                           TCounter;

        typedef typename Iterator<TText const, Standard>::Type      TTextIterator;
        typedef typename Iterator<TSA, Standard>::Type              TSAIterator;
        typedef typename Iterator<TCounter, Standard>::Type         TCntIterator;
        typedef typename Iterator<TCounter const, Standard>::Type   TConstCntIterator;
        typedef typename Size<TText>::Type                          TTextSize;

        TText const &text = indexText(index);
        TCounter const &tempSA = index.tempSA;
        TCounter &occ = index.tempOcc;
        TCounter &bound = index.tempBound;

        // 1. clear counters and copy SA to temporary SA
        arrayFill(begin(occ, Standard()), end(occ, Standard()), 0);
        index.tempSA = infix(indexSA(index), left, right);

        // 2. count characters
        index.sentinelBound = 0;
        index.sentinelOcc = _wotdCountChars(occ, text, tempSA, prefixLen);

        // 3. cumulative sum
        TSize requiredSize = 0;
        if (index.interSentinelNodes) {
            if (index.sentinelOcc != 0)
                requiredSize = (index.sentinelOcc > 1)? 2: 1;    // insert *one* $-edge for all $_i suffices
        } else
            requiredSize = index.sentinelOcc;                    // insert each $_i suffix one-by-one

        requiredSize += _wotdCummulativeSum(bound, occ, left + index.sentinelOcc);
        index.sentinelBound = left;

        // 4. fill suffix array
        {
            TSA &sa = indexSA(index);
            TSAIterator saBeg = begin(sa, Standard());
            TCntIterator boundBeg = begin(bound, Standard());

            TTextIterator itText = begin(text, Standard()) + prefixLen;
            TConstCntIterator itSA = begin(tempSA, Standard());
            TConstCntIterator itSAEnd = end(tempSA, Standard());
            TTextSize textLength = length(text) - prefixLen;
            for(; itSA != itSAEnd; ++itSA)
            {
                TTextSize saValue = *itSA;
                if (textLength > saValue)
                    *(saBeg + (*(boundBeg + ordValue(*(itText + saValue))))++) = saValue;
                else
                    if (textLength == saValue)
                        *(saBeg + index.sentinelBound++) = saValue;
            }
        }

        return requiredSize;
    }

    // multiple sequences
    template < typename TText, typename TSpec, typename TIndexSpec, typename TBeginPos, typename TEndPos, typename TSize >
    TSize _sortWotdBucket(
        Index<StringSet<TText, TSpec>, TIndexSpec> &index,
        TBeginPos left,
        TEndPos right,
        TSize prefixLen)
    {
        typedef Index<StringSet<TText, TSpec>, TIndexSpec>            TIndex;
        typedef typename Fibre<TIndex, WotdSA >::Type                TSA;
        typedef typename TIndex::TCounter                            TCounter;
        typedef typename TIndex::TTempSA                            TTempSA;

        typedef typename Iterator<TText const, Standard>::Type        TTextIterator;
        typedef typename Iterator<TSA, Standard>::Type                TSAIterator;
        typedef typename Iterator<TTempSA const, Standard>::Type    TTempSAIterator;
        typedef typename Iterator<TCounter, Standard>::Type            TCntIterator;
        typedef typename Size<TText>::Type                            TTextSize;

        StringSet<TText, TSpec> const &stringSet = indexText(index);
        TTempSA const &tempSA = index.tempSA;
        TCounter &occ = index.tempOcc;
        TCounter &bound = index.tempBound;

        // 1. clear counters and copy SA to temporary SA
        TCntIterator occBeg = begin(occ, Standard());

        arrayFill(occBeg, end(occ, Standard()), 0);
        index.tempSA = infix(indexSA(index), left, right);

        // 2. count characters
        index.sentinelBound = 0;
        index.sentinelOcc = _wotdCountChars(occ, stringSet, tempSA, prefixLen);

        // 3. cumulative sum
        TSize requiredSize = 0;
        if (index.interSentinelNodes) {
            if (index.sentinelOcc != 0)
                requiredSize = (index.sentinelOcc > 1)? 2: 1;    // insert *one* $-edge for all $_i suffices
        } else
            requiredSize = index.sentinelOcc;                    // insert each $_i suffix one-by-one

        requiredSize += _wotdCummulativeSum(bound, occ, left + index.sentinelOcc);
        index.sentinelBound = left;

        // 4. fill suffix array
        {
            if (empty(stringSet))
                return requiredSize;

            TSA &sa = indexSA(index);
            TSAIterator saBeg = begin(sa, Standard());
            TCntIterator boundBeg = begin(bound, Standard());

            TTextIterator itText = begin(front(stringSet), Standard());
            TTempSAIterator itSA = begin(tempSA, Standard());
            TTempSAIterator itSAEnd = end(tempSA, Standard());
            TTextSize textLength = 0;
            unsigned lastSeqSeen = (unsigned)-1;
            Pair<unsigned, TTextSize> lPos;
            for(; itSA != itSAEnd; ++itSA)
            {
                posLocalize(lPos, *itSA, stringSetLimits(index));
                if (lastSeqSeen != getSeqNo(lPos))
                {
                    lastSeqSeen = getSeqNo(lPos);

                    // shift textBegin and textLength by prefixLen
                    textLength = length(stringSet[lastSeqSeen]) - prefixLen;
                    itText = begin(stringSet[lastSeqSeen], Standard()) + prefixLen;
                }
                if (textLength > getSeqOffset(lPos))
                    *(saBeg + (*(boundBeg + ordValue(*(itText + getSeqOffset(lPos)))))++) = *itSA;
                else
                    if (textLength == getSeqOffset(lPos))
                        *(saBeg + index.sentinelBound++) = *itSA;
            }
        }

        return requiredSize;
    }





    template < typename TSA, typename TText >
    typename Size<TText>::Type
    _bucketLcp(TSA const &sa, TText const &text)
    {
        typedef typename Iterator<TText const, Standard>::Type    TTextIterator;
        typedef typename Iterator<TSA const, Standard>::Type    TSAIterator;
        typedef typename Value<TText>::Type                        TValue;
        typedef typename Size<TText>::Type                        TTextSize;

        TTextSize prefixLen = 0;

        if (length(sa) < 2) return prefixLen;

        TTextIterator itText = begin(text, Standard());
        TSAIterator itSAEnd = end(sa, Standard());
        TTextSize textLength = length(text);

        do {
            TSAIterator itSA = begin(sa, Standard());
            TTextSize sa = *itSA;
            if (textLength == sa) return prefixLen;

            TValue c = *(itText + sa);
            for(++itSA; itSA != itSAEnd; ++itSA) {
                sa = *itSA;
                if (textLength == sa || c != *(itText + sa))
                    return prefixLen;
            }
            ++prefixLen; --textLength;
            ++itText;
        } while (true);
    }

    template < typename TSA, typename TText, typename TSize >
    typename Size<TText>::Type
    _bucketLcp(TSA const &sa, TText const &text, TSize prefixLen)
    {
        typedef typename Iterator<TText const, Standard>::Type    TTextIterator;
        typedef typename Iterator<TSA const, Standard>::Type    TSAIterator;
        typedef typename Value<TText>::Type                        TValue;
        typedef typename Size<TText>::Type                        TTextSize;

        if (length(sa) < 2) return prefixLen;

        TTextIterator itText = begin(text, Standard()) + prefixLen;
        TSAIterator itSAEnd = end(sa, Standard());
        TTextSize textLength = length(text) - prefixLen;

        do {
            TSAIterator itSA = begin(sa, Standard());
            TTextSize sa = *itSA;
            if (textLength == sa) return prefixLen;

            TValue c = *(itText + sa);
            for(++itSA; itSA != itSAEnd; ++itSA) {
                sa = *itSA;
                if (textLength == sa || *(itText + sa) != c)
                    return prefixLen;
            }
            ++prefixLen; --textLength;
            ++itText;
        } while (true);
    }

    template < typename TSA, typename TText, typename TSpec, typename TSize >
    typename Size<TText>::Type
    _bucketLcp(TSA const &sa, StringSet<TText, TSpec> const &stringSet, TSize prefixLen)
    {
        typedef typename Iterator<TText const, Standard>::Type    TTextIterator;
        typedef typename Iterator<TSA const, Standard>::Type    TSAIterator;
        typedef typename Value<TText>::Type                        TValue;
        typedef typename Size<TText>::Type                        TTextSize;

        if (length(sa) < 2) return prefixLen;

        TSAIterator itSAEnd = end(sa, Standard());
        TTextIterator itText;
        TTextSize textLength;

        Pair<unsigned, TTextSize> lPos;
        do {
            TSAIterator itSA = begin(sa, Standard());
            posLocalize(lPos, *itSA, stringSetLimits(stringSet));

            unsigned lastSeqSeen = getSeqNo(*itSA);
            textLength = length(stringSet[lastSeqSeen]) - prefixLen;
            if (textLength == getSeqOffset(lPos)) return prefixLen;

            itText = begin(stringSet[lastSeqSeen], Standard()) + prefixLen;
            TValue c = *(itText + getSeqOffset(*itSA));
            for(++itSA; itSA != itSAEnd; ++itSA)
            {
                posLocalize(lPos, *itSA, stringSetLimits(stringSet));

                if (lastSeqSeen != getSeqNo(lPos))
                {
                    lastSeqSeen = getSeqNo(lPos);

                    // shift textBegin and textLength by prefixLen
                    textLength = length(stringSet[lastSeqSeen]) - prefixLen;
                    itText = begin(stringSet[lastSeqSeen], Standard()) + prefixLen;
                }

                if (textLength == getSeqOffset(lPos) || c != *(itText + getSeqOffset(lPos)))
                    return prefixLen;
            }
            ++prefixLen; --textLength;
            ++itText;
        } while (true);
    }


    template <typename TText, typename TSpec, typename TPos>
    inline TPos
    _getNodeLP(
        Index<TText, IndexWotd<TSpec> > const &index,
        TPos pos)
    {
        TPos w0 = dirAt(pos, index);
        if (w0 & index.LEAF)
            return w0 & index.BITMASK0;

        TPos w1 = dirAt(pos + 1, index);
        if (w1 & index.UNEVALUATED)
            return saAt(w0 & index.BITMASK0, index);
        else
            return w0 & index.BITMASK0;
    }

    // store buckets into directory
    // ORIGINAL VERSION: storing SA entries and topology links in Dir
    template <typename TText, typename TPos>
    inline void
    _storeWotdChildren(
        Index<TText, IndexWotd<WotdOriginal> > &index,
        TPos dirOfs)
    {
        typedef Index<TText, IndexWotd<WotdOriginal> >        TIndex;
        typedef typename Fibre<TIndex, WotdDir>::Type        TDir;
        typedef typename Iterator<TDir, Standard>::Type        TDirIterator;
        typedef typename Size<TDir>::Type                    TDirSize;
        typedef typename TIndex::TCounter                    TCounter;
        typedef typename Iterator<TCounter, Standard>::Type    TCntIterator;

        typedef typename Value<TCounter>::Type                TValue;

        TDirIterator itDir = begin(indexDir(index), Standard()) + dirOfs;
        TDirIterator itDirEnd = end(indexDir(index), Standard());
        TDirIterator itPrev = itDirEnd;

        TCntIterator it = begin(index.tempOcc, Standard());
        TCntIterator bit = begin(index.tempBound, Standard());
        TCntIterator itEnd = end(index.tempOcc, Standard());

        TValue occ;
        if (index.sentinelOcc != 0)
        {
            if (index.sentinelOcc > 1 && index.interSentinelNodes)    // occurs on multiseqs
            {
                itPrev = itDir;
                *itDir = index.sentinelBound - index.sentinelOcc;    ++itDir;
                *itDir = index.sentinelBound | index.UNEVALUATED;    ++itDir;
            } else
                for (TDirSize d = index.sentinelBound - index.sentinelOcc; d != index.sentinelBound; ++d)
                {
                    itPrev = itDir;
                    *itDir = saAt(d, index) | index.LEAF;            ++itDir;
                }
        }
        for (; it != itEnd; ++it, ++bit)
        {
            if ((occ = *it) == 0) continue;
            if (occ > 1) {
                itPrev = itDir;
                *itDir = *bit - occ;                             ++itDir;
                *itDir = *bit | index.UNEVALUATED;                ++itDir;
            } else {
                itPrev = itDir;
                *itDir = saAt(*bit - occ, index) | index.LEAF;    ++itDir;
            }
        }

        // mark the last child
        if (itPrev != itDirEnd)
            *itPrev |= index.LAST_CHILD;
    }

    // store buckets into directory
    // MODIFIED VERSION: storing SA links and topology links in Dir
    template <typename TText, typename TSpec, typename TSize>
    inline void
    _storeWotdChildren(
        Index<TText, IndexWotd<TSpec> > &index,
        TSize dirOfs,
        TSize lcp)
    {
        typedef Index<TText, IndexWotd<TSpec> >            TIndex;
        typedef typename Fibre<TIndex, WotdDir>::Type        TDir;
        typedef typename Iterator<TDir, Standard>::Type        TDirIterator;
        typedef typename Size<TDir>::Type                    TDirSize;
        typedef typename TIndex::TCounter                    TCounter;
        typedef typename Iterator<TCounter, Standard>::Type    TCntIterator;

        typedef typename Value<TCounter>::Type                TValue;

        TDirIterator itDirBegin = begin(indexDir(index), Standard()) + dirOfs;
        TDirIterator itDirEnd = end(indexDir(index), Standard());
        TDirIterator itDir = itDirBegin;
        TDirIterator itPrev = itDirEnd;

        TCntIterator it = begin(index.tempOcc, Standard());
        TCntIterator bit = begin(index.tempBound, Standard());
        TCntIterator itEnd = end(index.tempOcc, Standard());

        TValue occ;
        if (index.sentinelOcc != 0)
        {
            if (index.sentinelOcc > 1 && index.interSentinelNodes)    // occurs on multiseqs
            {
                itPrev = itDir;
                *itDir = index.sentinelBound - index.sentinelOcc;    ++itDir;
                *itDir = index.sentinelBound | index.UNEVALUATED;    ++itDir;
            } else
                for (TDirSize d = index.sentinelBound - index.sentinelOcc; d != index.sentinelBound; ++d)
                {
                    itPrev = itDir;
                    *itDir = d | index.LEAF;                        ++itDir;
                }
        }
        for (; it != itEnd; ++it, ++bit)
        {
            if ((occ = *it) == 0) continue;
            if (occ > 1) {
                itPrev = itDir;
                *itDir = *bit - occ;                 ++itDir;
                *itDir = *bit | index.UNEVALUATED;    ++itDir;
            } else {
                itPrev = itDir;
                *itDir = (*bit - occ) | index.LEAF;    ++itDir;
            }
        }

        // first child gets the mutual lcp value of the children (== parent repLength)
        *itDirBegin = ((*itDirBegin) & ~index.BITMASK0) | lcp;

        // mark the last child
        if (itPrev != itDirEnd)
            *itPrev |= index.LAST_CHILD;
    }


    template < typename TText, typename TSpec >
    inline typename Size< Index<TText, IndexWotd<WotdOriginal> > >::Type
    _wotdEvaluate(Iter< Index<TText, IndexWotd<WotdOriginal> >, VSTree<TSpec> > const &it)
    {
        typedef Index<TText, IndexWotd<WotdOriginal> >    TIndex;
        typedef typename Size<TIndex>::Type                TSize;

        TIndex &index = const_cast<TIndex&>(container(it));
        TSize pos = value(it).node;
        TSize w1 = dirAt(pos + 1, index);

        // test for evaluation
        if (w1 & index.UNEVALUATED)
        {
            TSize w0 = dirAt(pos, index);
            TSize lp = saAt(w0 & index.BITMASK0, index);
            TSize dst = length(indexDir(index));

            TSize size = _sortWotdBucket(
                index,
                w0 & index.BITMASK0,
                w1 & index.BITMASK1,
                parentEdgeLength(it));

            resize(indexDir(index), dst + size, Generous());
            _storeWotdChildren(index, dst);

            // mark nodes with solely empty child edges
            w1 = dst;
            if (index.sentinelOcc > 0)
            {
                TSize sentinelSize = index.sentinelOcc;
                if (index.interSentinelNodes && sentinelSize > 2)
                    sentinelSize = 2;
                if (size == sentinelSize) w1 |= index.SENTINELS;
            }

            assert(!(index.sentinelOcc == 1 && size == 1));

            dirAt(pos, index)     = (w0 & ~index.BITMASK0) | lp;
            dirAt(pos + 1, index) = w1;
        }

        return w1;
    }

    template < typename TText, typename TIndexSpec, typename TSpec >
    inline typename Size< Index<TText, IndexWotd<TIndexSpec> > >::Type
    _wotdEvaluate(Iter< Index<TText, IndexWotd<TIndexSpec> >, VSTree<TSpec> > const &it)
    {
        typedef Index<TText, IndexWotd<TIndexSpec> >    TIndex;
        typedef typename Size<TIndex>::Type                TSize;

        TIndex &index = const_cast<TIndex&>(container(it));
        TSize pos = value(it).node;
        TSize w1 = dirAt(pos + 1, index);

        // test for evaluation
        if (w1 & index.UNEVALUATED)
        {
            TSize dst = length(indexDir(index));

            TSize size = _sortWotdBucket(
                index,
                value(it).range.i1,
                w1 & index.BITMASK1,
                repLength(it));
/*
            if (globalDumpFlag)
            {
                std::cerr << '"' << representative(it) << '"' << std::endl;
                for (int i=0;i<length(getOccurrences(it));++i)
                    std::cerr << getOccurrences(it)[i]<<'\t'<<suffix(indexText(index),getOccurrences(it)[i])<<std::endl;
//                _dumpFreq(index);
            }
*/
            resize(indexDir(index), dst + size, Generous());
            _storeWotdChildren(index, dst, repLength(it));

            // mark nodes with solely empty child edges
            w1 = dst;
            if (index.sentinelOcc > 0)
            {
                TSize sentinelSize = index.sentinelOcc;
                if (index.interSentinelNodes && sentinelSize > 2)
                    sentinelSize = 2;
                if (size == sentinelSize) w1 |= index.SENTINELS;
            }

            dirAt(pos + 1, index) = w1;
        }

        return w1;
    }


    template <typename TText, typename TSpec>
    inline void
    _dump(Index<TText, IndexWotd<TSpec> > &index)
    {
        typedef Index<TText, IndexWotd<TSpec> >            TIndex;
        typedef typename Fibre<TIndex, WotdDir>::Type        TDir;
        typedef typename Value<TDir>::Type                    TDirValue;

        std::cout << "  Dir (wotd)" << std::endl;
        for(unsigned i=0; i < length(indexDir(index)); ++i) {
            TDirValue d = indexDir(index)[i];
            std::cout << i << ":  " << (d & index.BITMASK0);
            if (d & index.LEAF)            std::cout << "  (Leaf/Uneval)";
            if (d & index.LAST_CHILD)    std::cout << "  (LastChild/SENTINELS)";
            std::cout << std::endl;
        }

        std::cout << std::endl << "  SA" << std::endl;
        for(unsigned i=0; i < length(indexSA(index)); ++i)
            std::cout << i << ":  " << indexSA(index)[i] << "  " << suffix(indexText(index), indexSA(index)[i]) << std::endl;

        std::cout << std::endl;

    }

//////////////////////////////////////////////////////////////////////////////
// _goDownChar

    template < typename TText, class TSpec, typename TValue >
    inline bool _goDownChar(
        Iter<Index<TText, IndexWotd<WotdOriginal> >, VSTree< TopDown<TSpec> > > &it,
        TValue c)
    {
        typedef Index<TText, IndexWotd<TSpec> >    TIndex;
        typedef typename Value<TIndex>::Type        TIndexValue;

        bool sorted = false;
        if (!goDown(it)) return false;
        do {
            if (parentEdgeLength(it) != 0) {
                TIndexValue edgeChar = parentEdgeLabel(it)[0];
                if (edgeChar == c) return true;        // the edge is found
                if (sorted && edgeChar > c) break;    // too far (except the first one,
            }                                        // child edges are sorted)
            sorted = true;
        } while (goRight(it));
        _goUp(it);
        return false;
    }

    template < typename TText, class TIndexSpec, class TSpec, typename TValue >
    inline bool _goDownChar(
        Iter<Index<TText, IndexWotd<TIndexSpec> >, VSTree< TopDown<TSpec> > > &it,
        TValue c)
    {
        typedef Index<TText, IndexWotd<TSpec> >    TIndex;
        typedef typename Value<TIndex>::Type        TIndexValue;

        if (!goDown(it)) return false;
        do {
            if (parentEdgeLength(it) != 0) {
                TIndexValue edgeChar = parentEdgeLabel(it)[0];
                if (edgeChar == c) return true;    // the edge is found
                if (edgeChar > c) break;        // too far (child edges are sorted)
            }
        } while (goRight(it));
        _goUp(it);
        return false;
    }

/*
    template < typename TText, typename TSpec, typename TValue >
    inline bool
    _getNodeByChar(
        Iter< Index<TText, IndexWotd<TSpec> >, VSTree<TSpec> > const &it,
        TValue c,
        typename VertexDescriptor< Index<TText, IndexWotd<TSpec> > >::Type &childDesc)
    {
        typedef Index<TText, IndexWotd<TSpec> >            TIndex;
        typedef typename Fibre<TIndex, WotdDir>::Type        TDir;
        typedef typename Fibre<TIndex, WotdSA>::Type        TSA;

        typedef typename Value<TSA>::Type                    TSAValue;
        typedef typename Value<TDir>::Type                    TDirValue;

        typename Size<TIndex>::Type len = length(index);
        typename VertexDescriptor<TIndex>::Type    desc;

        TSAValue pos = _firstSuffixOfBucket(index, value(it).node);
        while (pos == len || value < (c = textAt(index, pos + value.i2))) {
            value.node += (dirAt(value.node, index) & index.LEAF)? 1: 2;
            pos = _firstSuffixOfBucket(index, value.node);
        }
        assert(pos != len);

        return c == value;
    }
*/

//////////////////////////////////////////////////////////////////////////////
// interface for automatic index creation

    template <typename TText, typename TSpec>
    inline void _wotdCreateFirstLevel(Index<TText, IndexWotd<TSpec> > &index)
    {
        typedef Index<TText, IndexWotd<TSpec> > TIndex;
        typedef typename Value<TIndex>::Type    TValue;
        typedef typename Size<TIndex>::Type     TSize;

        resize(index.tempOcc, ValueSize<TValue>::VALUE + 1, Exact());
        resize(index.tempBound, ValueSize<TValue>::VALUE + 1, Exact());

        TSize size;
        if (empty(indexSA(index)))
        {
            resize(indexSA(index), length(indexRawText(index)), Exact());
            size = _sortFirstWotdBucket(index);
        } else
        {
            size = _sortWotdBucket(index, 0, length(indexSA(index)), 0);
        }

        if (size > 0)
        {
            resize(indexDir(index), size + 2, Generous());
            _storeWotdChildren(index, 2, 0);

            // mark nodes with solely empty child edges
            TSize w1 = 2;
            if (index.sentinelOcc > 0)
            {
                TSize sentinelSize = index.sentinelOcc;
                if (index.interSentinelNodes && sentinelSize > 2)
                    sentinelSize = 2;
                if (size == sentinelSize) w1 |= index.SENTINELS;
            }


            dirAt(0, index) = 0 | index.LAST_CHILD;
            dirAt(1, index) = w1;

        } else {
            resize(indexDir(index), 1);
            dirAt(0, index) = 0 | index.LAST_CHILD | index.LEAF;
        }
    }

    template <typename TText, typename TSpec>
    inline bool indexCreate(Index<TText, IndexWotd<TSpec> > &index, WotdDir const, Default const)
    {
        _wotdCreateFirstLevel(index);
        return true;
    }

//////////////////////////////////////////////////////////////////////////////
// clear

    template < typename TText, typename TSpec >
    inline void clear(Index<TText, IndexWotd<TSpec> > &index)
    {
        clear(getFibre(index, WotdSA()));
        clear(getFibre(index, WotdDir()));
    }

//////////////////////////////////////////////////////////////////////////////
// open

    template < typename TText, typename TSpec >
    inline bool open(
        Index< TText, IndexWotd<TSpec> > &index,
        const char *fileName,
        int openMode)
    {
        typedef Index<TText, IndexWotd<TSpec> > TIndex;
        typedef typename Value<TIndex>::Type    TValue;

        String<char> name;

        name = fileName;    append(name, ".txt");
        if ((!open(getFibre(index, WotdText()), toCString(name), openMode)) &&
            (!open(getFibre(index, WotdText()), fileName, openMode))) return false;

        name = fileName;    append(name, ".sa");
        if (!open(getFibre(index, WotdSA()), toCString(name), openMode)) return false;
        name = fileName;    append(name, ".dir");

        if (!open(getFibre(index, WotdDir()), toCString(name), openMode)) return false;

        if (!empty(getFibre(index, WotdDir())))
        {
            resize(index.tempOcc, ValueSize<TValue>::VALUE + 1);
            resize(index.tempBound, ValueSize<TValue>::VALUE + 1);
        }

        return true;
    }
    template < typename TText, typename TSpec >
    inline bool open(
        Index< TText, IndexWotd<TSpec> > &index,
        const char *fileName)
    {
        return open(index, fileName, OPEN_RDONLY);
    }


//////////////////////////////////////////////////////////////////////////////
// save

    template < typename TText, typename TSpec >
    inline bool save(
        Index< TText, IndexWotd<TSpec> > &index,
        const char *fileName,
        int openMode)
    {
        String<char> name;

        name = fileName;    append(name, ".txt");
        if ((!save(getFibre(index, WotdText()), toCString(name), openMode)) &&
            (!save(getFibre(index, WotdText()), fileName, openMode))) return false;

        name = fileName;    append(name, ".sa");
        if (!save(getFibre(index, WotdSA()), toCString(name), openMode)) return false;
        name = fileName;    append(name, ".dir");

        if (!save(getFibre(index, WotdDir()), toCString(name), openMode)) return false;

        return true;
    }
    template < typename TText, typename TSpec >
    inline bool save(
        Index< TText, IndexWotd<TSpec> > &index,
        const char *fileName)
    {
        return save(index, fileName, OPEN_WRONLY | OPEN_CREATE);
    }
}

#endif //#ifndef SEQAN_HEADER_...
