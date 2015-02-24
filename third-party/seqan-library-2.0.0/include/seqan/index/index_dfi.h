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

#ifndef SEQAN_HEADER_INDEX_DFI_H
#define SEQAN_HEADER_INDEX_DFI_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
// simple struct to store frequency vectors

    struct DfiEntry_
    {
        unsigned                        lastSeqSeen;
//        String<unsigned, Array<2> >        freq;
        String<unsigned>                freq;
    };


//////////////////////////////////////////////////////////////////////////////
// constant frequency predicate

    template <bool Result>
    struct DfiPredDefault_
    {
        inline bool operator()(DfiEntry_ const &) const {
            return Result;
        }
    };


//////////////////////////////////////////////////////////////////////////////
// Dfi - The Deferred Frequency Index

/*!
 * @class IndexDfi
 * @extends IndexWotd
 * @headerfile <seqan/index.h>
 *
 * @brief The Deferred Frequency Index (see Weese and Schulz, "Efficient string mining under constraints via
 *        the deferred frequency index").
 *
 * @signature template <typename TText, typename TPredHull, typename TPred>
 *            class Index<TText, IndexWotd< Dfi<TPredHull, TPred> > >;
 *
 * @tparam TText     The @link TextConcept text @endlink.
 * @tparam TPred     An arbitrary frequency predicate.
 * @tparam TPredHull A monotonic hull of <tt>TPred</tt>
 *
 * This index is based on a lazy suffix tree (see @link IndexWotd @endlink).  All <tt>TPredHull</tt> sufficing
 * nodes can be iterated using a @link TopDownIterator @endlink.  To iterate the exact solution set of <tt>TPred</tt>,
 * use a @link TopDownHistoryIterator @endlink of this index.
 */

/*!
 * @defgroup DfiIndexFibres Dfi Index Fibres
 * @brief Tag to select a specific fibre (e.g. table, object, ...) of an @link
 *        IndexDfi @endlink index.
 *
 * These tags can be used to get @link Fibre @endlink of an @link IndexDfi @endlink.
 *
 * @see Fibre
 * @see Index#getFibre
 * @see IndexDfi
 *
 * @tag DfiIndexFibres#DfiDir
 * @brief The child table.
 *
 * @tag DfiIndexFibres#DfiRawSA
 * @brief The raw suffix array.
 *
 * @tag DfiIndexFibres#DfiText
 * @brief The original text the index should be based on.
 *
 * @tag DfiIndexFibres#DfiRawText
 * @brief The raw text the index is really based on.
 *
 * @tag DfiIndexFibres#DfiSA
 * @brief The suffix array.
 */

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexDfi#indexSA
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>getFibre(.., DfiSA)</tt>.
 *
 * @signature TSa indexSA(index);
 *
 * @param[in] index The @link IndexDfi @endlink object holding the fibre.
 *
 * @return TSa A reference to the @link DfiIndexFibres#DfiSA @endlink fibre (partially sorted suffix array).
 */

/*!
 * @fn IndexDfi#indexDir
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>getFibre(.., DfiDir())</tt>.
 * @signature TFibre indexDir(index);
 *
 * @param[in] index The @link IndexDfi @endlink object holding the fibre.
 *
 * @return TFibre A reference to the @link DfiIndexFibres#DfiDir @endlink fibre (tree structure).
 */

/*!
 * @fn IndexDfi#saAt
 * @headerfile <seqan/index.h>
 * @note Advanced functionality, not commonly used.
 * @brief Shortcut for <tt>value(indexSA(..), ..)</tt>.
 *
 * @signature TValue saAt(position, index);
 *
 * @param[in] index The @link IndexDfi @endlink object holding the fibre.
 * @param[in] position A position in the array on which the value should be accessed.
 *
 * @return TValue A reference or proxy to the value in the @link DfiIndexFibres#DfiSA @endlink fibre.
 *                To be more precise, a reference to a position containing a value of type
 *                @link SAValue @endlink is returned (or a proxy).
 */

/*!
 * @fn IndexDfi#dirAt
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>value(indexDir(index), position)</tt>.
 *
 * @signature TFibre dirAt(position, index);
 *
 * @param[in] index    The @link IndexDfi @endlink object holding the fibre.
 * @param[in] position A position in the array on which the value should be accessed.
 *
 * @return TFibre A reference to the @link DfiIndexFibres#DfiDir @endlink fibre.
 */

        template <
        typename TPredHull = DfiPredDefault_<true>,
        typename TPred = DfiPredDefault_<true>
    >
    struct Dfi;

    template <
        typename TObject,
        typename TPredHull,
        typename TPred
    >
    class Index<TObject, IndexWotd< Dfi<TPredHull, TPred> > >:
        public Index<TObject, IndexWotd<> >
    {
    public:

        typedef Index<TObject, IndexWotd<> >    TBase;

        // extending base class
        typedef typename TBase::TText    TText;
        typedef typename TBase::TValue    TValue;
        typedef typename TBase::TSize    TSize;

        using TBase::LEAF;
        using TBase::LAST_CHILD;
        using TBase::UNEVALUATED;
        using TBase::SENTINELS;

        // frequency strings for Dfi
        typedef DfiEntry_                        TDFIEntry;
        typedef String<
            TDFIEntry,
            Array<ValueSize<TValue>::VALUE> >    TDFIEntries;
        typedef String<unsigned>                TDFIDatasets;

        // 1st word flags
        static TSize const DFI_PRED_HULL    = (TSize)1 << (BitsPerValue<TSize>::VALUE - 3); // this node fulfills all monotonic frequency predicates (e.g. min_freq)
        static TSize const DFI_PRED            = (TSize)1 << (BitsPerValue<TSize>::VALUE - 4); // this node fulfills all frequency predicates (e.g. emerging, minmax_freq)
        static TSize const DFI_PARENT_FREQ    = (TSize)1 << (BitsPerValue<TSize>::VALUE - 5); // this node has the same frequency as its parent

        static TSize const BITMASK0            = ~(LEAF | LAST_CHILD | DFI_PRED_HULL | DFI_PRED | DFI_PARENT_FREQ);
        static TSize const BITMASK1            = ~(UNEVALUATED | SENTINELS);

        TDFIEntry        nodeEntry;            // current nodes frequencies
        TDFIEntries        childEntry;            // child frequencies for each first letter of outgoing edges
        TDFIDatasets    ds;

        TPredHull        predHull;
        TPred            pred;

        Index() {}

        Index(Index &other):
            TBase((TBase &)other),
            childEntry(other.childEntry),
            ds(other.ds),
            predHull(other.predHull),
            pred(other.pred) {}

        Index(Index const &other):
            TBase((TBase const &)other),
            childEntry(other.childEntry),
            ds(other.ds),
            predHull(other.predHull),
            pred(other.pred) {}

        template <typename TText_>
        Index(TText_ &_text):
            TBase(_text) {}

        template <typename TText_>
        Index(TText_ &_text, TPredHull const &_predHull):
            TBase(_text),
            predHull(_predHull) {}

        template <typename TText_>
        Index(TText_ &_text, TPredHull const &_predHull, TPred const &_pred):
            TBase(_text),
            predHull(_predHull),
            pred(_pred) {}
    };


    template <
        typename TText,
        typename TPredHull,
        typename TPred,
        typename TSpec
    >
    inline bool nodePredicate(
        Iter<Index<TText, IndexWotd< Dfi<TPredHull, TPred> > >, TSpec> &it)
    {
        typedef Index<TText, IndexWotd< Dfi<TPredHull, TPred> > > TIndex;
        return (dirAt(value(it).node, container(it)) & TIndex::DFI_PRED) != 0;
    }

    template <
        typename TText,
        typename TPredHull,
        typename TPred,
        typename TSpec
    >
    inline bool nodeHullPredicate(
        Iter<Index<TText, IndexWotd< Dfi<TPredHull, TPred> > >, TSpec> &it)
    {
        typedef Index<TText, IndexWotd< Dfi<TPredHull, TPred> > > TIndex;
        return (dirAt(value(it).node, container(it)) & TIndex::DFI_PRED_HULL) != 0;
    }


//////////////////////////////////////////////////////////////////////////////
// we modify counting sort used in the wotd-index
// to count the frequencies in passing

    template < typename TText, typename TSpec, typename TPredHull, typename TPred >
    typename Size< Index<StringSet<TText, TSpec>, IndexWotd<Dfi<TPredHull, TPred> > > >::Type
    _sortFirstWotdBucket(Index<StringSet<TText, TSpec>, IndexWotd<Dfi<TPredHull, TPred> > > &index)
    {
    SEQAN_CHECKPOINT
        typedef Index<StringSet<TText, TSpec>, IndexWotd<Dfi<TPredHull, TPred> > >    TIndex;
        typedef typename Fibre<TIndex, WotdSA >::Type            TSA;
        typedef typename TIndex::TCounter                        TCounter;

        typedef typename TIndex::TDFIEntry                        TDFIEntry;
        typedef typename TIndex::TDFIDatasets                    TDFIDatasets;
        typedef typename Iterator<TDFIDatasets, Standard>::Type    TDFIDatasetsIterator;

        typedef typename Iterator<TText const, Standard>::Type    TTextIterator;
        typedef typename Iterator<TSA, Standard>::Type            TSAIterator;
        typedef typename Iterator<TCounter, Standard>::Type        TCntIterator;
        typedef typename Size<TText>::Type                        TSize;
        typedef typename Value<TText>::Type                        TValue;

        StringSet<TText, TSpec> const &stringSet = indexText(index);
        TCounter &occ = index.tempOcc;
        TCounter &bound = index.tempBound;

        // 1. clear counters and copy SA to temporary SA
        arrayFill(begin(occ, Standard()), end(occ, Standard()), 0);

        index.nodeEntry.lastSeqSeen = -1;
        for(unsigned j = 0; j < length(index.nodeEntry.freq); ++j)
            index.nodeEntry.freq[j] = 0;
        for(unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i)
        {
            TDFIEntry &childEntry = index.childEntry[i];
            childEntry.lastSeqSeen = -1;
            for(unsigned j = 0; j < length(childEntry.freq); ++j)
                childEntry.freq[j] = 0;
        }

        // 2. count characters
        _wotdCountChars(occ, stringSet);

        // 3. cummulative sum
        TSize requiredSize = _wotdCummulativeSum(bound, occ);

        // 4. fill suffix array
        unsigned dsNo = 0;
        TDFIDatasetsIterator currentDS = begin(index.ds, Standard()) + 1;
        for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo)
        {
            // search for surrounding dataset
            while (seqNo >= *currentDS) {
                ++dsNo;
                ++currentDS;
            }

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
                unsigned ord = ordValue(*itText);
                TDFIEntry &childEntry = index.childEntry[ord];
                // new sequence is seen for <ord> character
                // -> increment frequency of current dataset
                if (childEntry.lastSeqSeen != seqNo)
                {
                    childEntry.lastSeqSeen = seqNo;
                    ++childEntry.freq[dsNo];
                }

                *(saBeg + (*(boundBeg + ord))++) = localPos;
                assignValueI2(localPos, getValueI2(localPos) + 1);
            }
        }
        index.sentinelOcc = 0;
        index.sentinelBound = 0;

        return requiredSize;
    }

    // sort bucket using radixsort
    // - all buckets are in lexicographical order
    // - SA[left,right) contains real SA entries (the beginning positions of the suffices)
    template < typename TText, typename TSpec, typename TPredHull, typename TPred, typename TBeginPos, typename TEndPos, typename TSize >
    TSize _sortWotdBucket(
        Index<StringSet<TText, TSpec>, IndexWotd<Dfi<TPredHull, TPred> > > &index,
        TBeginPos left,
        TEndPos right,
        TSize prefixLen)
    {
    SEQAN_CHECKPOINT
        typedef Index<StringSet<TText, TSpec>, IndexWotd<Dfi<TPredHull, TPred> > >    TIndex;
        typedef typename Fibre<TIndex, WotdSA >::Type                TSA;
        typedef typename TIndex::TCounter                            TCounter;
        typedef typename TIndex::TTempSA                            TTempSA;
        typedef typename TIndex::TDFIEntry                            TDFIEntry;
        typedef typename TIndex::TDFIDatasets                        TDFIDatasets;
        typedef typename Iterator<TDFIDatasets, Standard>::Type        TDFIDatasetsIterator;

        typedef typename Iterator<TText const, Standard>::Type        TTextIterator;
        typedef typename Iterator<TSA, Standard>::Type                TSAIterator;
        typedef typename Iterator<TTempSA const, Standard>::Type    TTempSAIterator;
        typedef typename Iterator<TCounter, Standard>::Type            TCntIterator;
        typedef typename Size<TText>::Type                            TTextSize;
        typedef typename Value<TText>::Type                            TValue;

        StringSet<TText, TSpec> const &stringSet = indexText(index);
        TTempSA const &tempSA = index.tempSA;
        TCounter &occ = index.tempOcc;
        TCounter &bound = index.tempBound;

        // 1. clear counters and copy SA to temporary SA
        TCntIterator occBeg = begin(occ, Standard());

        arrayFill(occBeg, end(occ, Standard()), 0);
        index.tempSA = infix(indexSA(index), left, right);

        index.nodeEntry.lastSeqSeen = -1;
        for(unsigned j = 0; j < length(index.nodeEntry.freq); ++j)
            index.nodeEntry.freq[j] = 0;
        for(unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i) {
            TDFIEntry &childEntry = index.childEntry[i];
            childEntry.lastSeqSeen = -1;
            for(unsigned  j = 0; j < length(childEntry.freq); ++j)
                childEntry.freq[j] = 0;
        }

        index.sentinelOcc = 0;
        index.sentinelBound = 0;

        // 2. count characters
        {
            TDFIDatasetsIterator currentDS = begin(index.ds, Standard()) + 1;
            TTextIterator itText = TTextIterator();
            TTempSAIterator itSA = begin(tempSA, Standard());
            TTempSAIterator itSAEnd = end(tempSA, Standard());
            TTextSize textLength = 0;
            unsigned lastSeqSeen = (unsigned)-1;
            unsigned dsNo = 0;
            Pair<unsigned, TTextSize> lPos;
            for (; itSA != itSAEnd; ++itSA)
            {
                posLocalize(lPos, *itSA, stringSetLimits(index));
                if (lastSeqSeen != getSeqNo(lPos))
                {
                    lastSeqSeen = getSeqNo(lPos);

                    // search for surrounding dataset
                    while (lastSeqSeen >= *currentDS)
                    {
                        ++dsNo;
                        ++currentDS;
                    }
                    ++index.nodeEntry.freq[dsNo];

                    // shift textBegin and textLength by prefixLen
                    textLength = length(stringSet[lastSeqSeen]) - prefixLen;
                    itText = begin(stringSet[lastSeqSeen], Standard()) + prefixLen;
                }
                if (textLength > getSeqOffset(lPos))
                {
                    unsigned ord = ordValue(*(itText + getSeqOffset(lPos)));
                    TDFIEntry &childEntry = index.childEntry[ord];
                    // new sequence is seen for <ord> character
                    // -> increment frequency of current dataset
                    if (childEntry.lastSeqSeen != lastSeqSeen)
                    {
                        childEntry.lastSeqSeen = lastSeqSeen;
                        ++childEntry.freq[dsNo];
                    }
                    ++*(occBeg + ord);
                } else
                    if (textLength == getSeqOffset(lPos)) ++index.sentinelOcc;
            }
        }

        // 3. cumulative sum
        TSize requiredSize = 0;
        if (index.interSentinelNodes) {
            if (index.sentinelOcc != 0)
                requiredSize = (index.sentinelOcc > 1)? 2: 1;    // insert *one* $-edge for all $_i suffices
        } else
            requiredSize = index.sentinelOcc;                    // insert each $_i suffix one-by-one

        requiredSize += _wotdCummulativeSum(bound, occ, left + index.sentinelOcc);
        index.sentinelBound = left;
/*
        std::cout << "$=" << index.sentinelOcc<<"@"<<index.sentinelBound << "\t";
        for(int i=0; i<length(occ);++i)
            if (occ[i])
                std::cout << i << "=" << occ[i]<<"@"<<bound[i] << "\t";
*/
        // 4. fill suffix array
        {
            TSA &sa = indexSA(index);
            TSAIterator saBeg = begin(sa, Standard());
            TCntIterator boundBeg = begin(bound, Standard());

            TTextIterator itText = TTextIterator();
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


    // store buckets into directory
    // storing SA links and topology links in Dir
    template <typename TText, typename TPredHull, typename TPred, typename TSize>
    inline void
    _storeWotdChildren(
        Index<TText, IndexWotd<Dfi<TPredHull, TPred> > > &index,
        TSize dirOfs,
        TSize lcp)
    {
    SEQAN_CHECKPOINT
        typedef Index<TText, IndexWotd<Dfi<TPredHull, TPred> > >    TIndex;

        typedef typename Fibre<TIndex, WotdDir>::Type        TDir;
        typedef typename TIndex::TCounter                    TCounter;
        typedef typename TIndex::TDFIEntries                TEntries;

        typedef typename Iterator<TDir, Standard>::Type        TDirIterator;
        typedef typename Size<TDir>::Type                    TDirSize;
        typedef typename Iterator<TCounter, Standard>::Type    TCntIterator;
        typedef typename Iterator<TEntries, Standard>::Type    TEntriesIterator;

        typedef typename Value<TCounter>::Type                TCntValue;
        typedef typename Value<TDir>::Type                    TDirValue;

        TDirIterator itDirBegin = begin(indexDir(index), Standard()) + dirOfs;
        TDirIterator itDirEnd = end(indexDir(index), Standard());
        TDirIterator itDir = itDirBegin;
        TDirIterator itPrev = itDirEnd;

        TCntIterator it = begin(index.tempOcc, Standard());
        TCntIterator bit = begin(index.tempBound, Standard());
        TCntIterator itEnd = end(index.tempOcc, Standard());
        TEntriesIterator itEntry = begin(index.childEntry, Standard());

        TCntValue occ;
        if (index.sentinelOcc != 0)
        {
            TDirValue orMask = (index.predHull(*itEntry))? index.DFI_PRED_HULL: 0;
            if (index.pred(*itEntry)) orMask |= index.DFI_PRED;

            if (index.sentinelOcc > 1 && index.interSentinelNodes)    // occurs on multiseqs
            {
                itPrev = itDir;
                *itDir = (index.sentinelBound - index.sentinelOcc) | orMask;    ++itDir;
                *itDir = index.sentinelBound | index.UNEVALUATED;                ++itDir;
            } else
                orMask |= index.LEAF;
                for (TDirSize d = index.sentinelBound - index.sentinelOcc; d != index.sentinelBound; ++d)
                {
                    itPrev = itDir;
                    *itDir = d | orMask;                                        ++itDir;
                }
        }

        for (; it != itEnd; ++it, ++bit, ++itEntry)
        {
            if ((occ = *it) == 0) continue;

            TDirValue orMask = (index.predHull(*itEntry))? index.DFI_PRED_HULL: 0;
            if (index.pred(*itEntry)) orMask |= index.DFI_PRED;
            if ((*itEntry).freq == index.nodeEntry.freq) orMask |= index.DFI_PARENT_FREQ;

            if (occ > 1) {
                itPrev = itDir;
                *itDir = (*bit - occ) | orMask;                    ++itDir;
                *itDir = *bit | index.UNEVALUATED;                ++itDir;
            } else {
                itPrev = itDir;
                *itDir = (*bit - occ) | index.LEAF | orMask;    ++itDir;
            }
        }

        // first child gets the mutual lcp value of the children (== parent repLength)
        *itDirBegin = ((*itDirBegin) & ~index.BITMASK0) | lcp;

        // mark the last child
        if (itPrev != itDirEnd)
            *itPrev |= index.LAST_CHILD;
    }

//////////////////////////////////////////////////////////////////////////////
// debug output

    template <typename TText, typename TPredHull, typename TPred>
    inline void
    _dump(Index<TText, IndexWotd< Dfi<TPredHull, TPred> > > &index)
    {
        typedef IndexWotd< Dfi<TPredHull, TPred> >        TSpec;
        typedef Index<TText, TSpec>                            TIndex;
        typedef typename Fibre<TIndex, WotdDir>::Type        TDir;
        typedef typename Value<TDir>::Type                    TDirValue;

        std::cout << "  Dir (wotd/Dfi)" << std::endl;
        for(unsigned i=0; i < length(indexDir(index)); ++i) {
            TDirValue d = indexDir(index)[i];
            std::cout << i << ":  " << (d & index.BITMASK0);
            if (d & index.LEAF)                std::cout << "  (Leaf/Uneval)";
            if (d & index.LAST_CHILD)        std::cout << "  (LastChild/SENTINELS)";
            if (d & index.DFI_PRED_HULL)    std::cout << "  (PRED_HULL)";
            if (d & index.DFI_PRED)            std::cout << "  (PRED)";
            if (d & index.DFI_PARENT_FREQ)    std::cout << "  (PARENT_FREQ)";
            std::cout << std::endl;
        }

        std::cout << std::endl << "  SA" << std::endl;
        for(unsigned i=0; i < length(indexSA(index)); ++i)
            std::cout << i << ":  " << indexSA(index)[i] << "  " << suffix(indexText(index), indexSA(index)[i]) << std::endl;

        std::cout << std::endl;
    }

    template <typename TText, typename TPredHull, typename TPred>
    inline void
    _dumpFreq(Index<TText, IndexWotd< Dfi<TPredHull, TPred> > > &index)
    {
        typedef Dfi<TPredHull, TPred> TSpec;
        typedef Index<TText, IndexWotd<TSpec> >                TIndex;
        typedef typename Value<TIndex>::Type                    TValue;

        std::cout << "  ParentF = (";
        for(unsigned d=0; d<length(index.nodeEntry.freq); ++d) {
            if (d>0) std::cout << ",";
            std::cout << index.nodeEntry.freq[d];
        }
        std::cout << ")" << std::endl;

        for(unsigned i=0; i<length(index.tempOcc); ++i)
            if (index.tempOcc[i] != 0) {
                std::cout << "  Freq[" << (TValue)i << "] = (";
                for(unsigned d=0; d<length(index.childEntry[i].freq); ++d) {
                    if (d>0) std::cout << ",";
                    std::cout << index.childEntry[i].freq[d];
                }
                std::cout << ")" << std::endl;
            }
    }

//////////////////////////////////////////////////////////////////////////////
// interface for automatic index creation

    template <typename TText, typename TPredHull, typename TPred>
    inline bool indexCreate(Index<TText, IndexWotd<Dfi<TPredHull, TPred> > > &index, WotdDir const, Default const)
    {
        typedef Index<TText, IndexWotd<Dfi<TPredHull, TPred> > >    TIndex;
        typedef typename Value<TIndex>::Type                            TValue;

        resize(index.childEntry, (unsigned) ValueSize<TValue>::VALUE);
        if (empty(index.ds)) {
            resize(index.ds, 2);
            index.ds[0] = 0;
            index.ds[1] = countSequences(index);
        }
        resize(index.nodeEntry.freq, length(index.ds) - 1, Exact());
        for(unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i)
            resize(index.childEntry[i].freq, length(index.ds) - 1, Exact());

        _wotdCreateFirstLevel(index);
        return true;
    }
}

#endif //#ifndef SEQAN_HEADER_...
