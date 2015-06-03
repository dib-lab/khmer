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
// Author: Tobias Rausch <rausch@embl.de>
// ==========================================================================

#ifndef SEQAN_CONSENSUS_CONSENSUS_REALIGN_H_
#define SEQAN_CONSENSUS_CONSENSUS_REALIGN_H_

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////////

template <typename TPos, typename TGapAnchor, typename TSpec, typename TGapPos>
inline void
removeGap(AlignedReadStoreElement<TPos, TGapAnchor, TSpec>& alignedRead,
          TGapPos const gapPos)
{
    typedef String<TGapAnchor> TGaps;
    typedef typename Iterator<TGaps, Standard>::Type TGapIter;
    if (gapPos < (TGapPos) alignedRead.beginPos) {
        --alignedRead.beginPos; --alignedRead.endPos;
    } else if (gapPos < (TGapPos) alignedRead.endPos) {
        --alignedRead.endPos;
        TGapIter gapIt = upperBoundGapAnchor(alignedRead.gaps, gapPos - alignedRead.beginPos, SortGapPos() );
        TGapIter gapItEnd = end(alignedRead.gaps, Standard());
        // Note: We might create empty gaps here
        for(;gapIt != gapItEnd; ++gapIt)
            --(gapIt->gapPos);
    }
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlignedReads, typename TSpec, typename TGapPos>
inline void
removeGap(String<TAlignedReads, TSpec>& alignedReadStore,
          TGapPos const gapPos)
{
    typedef String<TAlignedReads, TSpec> TAlignedReadStore;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type TAlignIter;

    TAlignIter alignIt = begin(alignedReadStore, Standard());
    TAlignIter alignItEnd = end(alignedReadStore, Standard());
    for(;alignIt != alignItEnd; ++alignIt)
        removeGap(*alignIt, gapPos);
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TPos, typename TGapAnchor, typename TSpec, typename TGapPos>
inline int
insertGap(AlignedReadStoreElement<TPos, TGapAnchor, TSpec>& alignedRead,
          TGapPos const gapPos)
{
    typedef String<TGapAnchor> TGaps;
    typedef typename Iterator<TGaps, Standard>::Type TGapIter;

    if (gapPos <= (TGapPos)alignedRead.beginPos) {
        ++alignedRead.beginPos; ++alignedRead.endPos;
        return 0;
    } else if (gapPos < (TGapPos)alignedRead.endPos) {
        ++alignedRead.endPos;
        TGapIter gapIt = lowerBoundGapAnchor(alignedRead.gaps, gapPos - alignedRead.beginPos, SortGapPos() );
        TGapIter gapItEnd = end(alignedRead.gaps, Standard());
        TGapPos insertPos = (gapPos - alignedRead.beginPos);
        if (gapIt == gapItEnd) {
            int gapLen = 0;
            if (gapItEnd != begin(alignedRead.gaps)) {
                --gapItEnd;
                gapLen = (int) gapItEnd->gapPos - (int) gapItEnd->seqPos;
            }
            appendValue(alignedRead.gaps, TGapAnchor(insertPos - gapLen, insertPos + 1), Generous());
        }
        else {
            int gapPrev = 0;
            if (gapIt != begin(alignedRead.gaps)) {
                TGapIter gapPrevious = gapIt;
                --gapPrevious;
                gapPrev = (int) gapPrevious->gapPos - (int) gapPrevious->seqPos;
            }
            // If gap is within an existing gap, extend this gap
            if (((TGapPos)(gapIt->gapPos - (((int) gapIt->gapPos - (int) gapIt->seqPos) - gapPrev)) <= insertPos) &&
                ((TGapPos)gapIt->gapPos >= insertPos)) {
                for(;gapIt != gapItEnd; ++gapIt)
                    ++(gapIt->gapPos);
            } else {
                // Otherwise, create a new gap
                TGapAnchor tmp = value(gapIt);
                ++tmp.gapPos;
                gapIt->gapPos = insertPos + 1;
                gapIt->seqPos = insertPos - gapPrev;
                do {
                    ++gapIt;
                    TGapAnchor newTmp;
                    if (gapIt != gapItEnd) {
                        newTmp = *gapIt;
                        ++newTmp.gapPos;
                        *gapIt = tmp;
                    } else appendValue(alignedRead.gaps, tmp, Generous() );
                    tmp = newTmp;
                } while (gapIt != gapItEnd);
            }
        }
        return 1;
    }
    return 0;
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlignedReads, typename TSpec, typename TGapPos>
inline int
insertGap(String<TAlignedReads, TSpec>& alignedReadStore,
          TGapPos const gapPos)
{
    typedef String<TAlignedReads, TSpec> TAlignedReadStore;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type TAlignIter;

    int numGaps = 0;
    TAlignIter alignIt = begin(alignedReadStore, Standard());
    TAlignIter alignItEnd = end(alignedReadStore, Standard());
    for(;alignIt != alignItEnd; ++alignIt)
        numGaps += insertGap(*alignIt, gapPos);
    return numGaps;
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TConsensus>
inline int
scoreConsensus(TConsensus& consensus)
{
    typedef typename Size<TConsensus>::Type TSize;
    typedef typename Iterator<TConsensus, Standard>::Type TConsIter;

    // Compute the score
    int score = 0;
    TConsIter itCons = begin(consensus, Standard() );
    TConsIter itConsEnd = end(consensus, Standard() );
    TSize maxCount = 0;
    TSize sumCount = 0;
    TSize tmp;
    for(;itCons != itConsEnd; ++itCons) {
        maxCount = 0; sumCount = 0;
        for(TSize i = 0; i < ValueSize<typename Value<TConsensus>::Type>::VALUE; ++i) {
            if ((tmp = (*itCons).count[i]) > maxCount) maxCount = tmp;
            sumCount += tmp;
        }
        score += (sumCount - maxCount);
    }
    return score;
}

//////////////////////////////////////////////////////////////////////////////////

// Perform one realignment round.
// TODO(holtgrew): Rename to reflect this more clearly.
// TODO(holtgrew): TConsensus/consensus are profiles, really.
template<typename TFragSpec, typename TConfig, typename TAlignedRead, typename TSpec, typename TConsensus, typename TScore, typename TMethod, typename TBandwidth>
void
reAlign(FragmentStore<TFragSpec, TConfig>& fragStore,
        String<TAlignedRead, TSpec>& contigReads,
        TConsensus& consensus,
        TScore& consScore,
        TMethod const rmethod,
        TBandwidth const bandwidth,
        bool includeReference,
        double & timeBeforeAlign,
        double & timeAlign,
        double & timeAfterAlign
        )
{
    typedef FragmentStore<TSpec, TConfig> TFragmentStore;
    typedef String<TAlignedRead, TSpec> TAlignedReadStore;
    typedef typename TFragmentStore::TReadPos TReadPos;
    typedef typename TFragmentStore::TReadSeq TReadSeq;
    typedef typename Size<TFragmentStore>::Type TSize;
    typedef typename TFragmentStore::TReadGapAnchor TGapAnchor;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type TAlignedReadIter;
    typedef typename Iterator<TConsensus, Standard>::Type TConsIter;

    // Initialization
    typedef typename Value<TConsensus>::Type TProfileChar;
    TSize gapPos = ValueSize<TProfileChar>::VALUE - 1;
//    int refId = length(fragStore.readSeqStore) - 1;

    // Remove each fragment and realign it to the profile.
//    TAlignedReadIter beg = begin(contigReads, Standard());
    TAlignedReadIter alignIt = begin(contigReads, Standard());
    TAlignedReadIter alignItEnd = end(contigReads, Standard());
    if (includeReference)
        --alignItEnd;
    TConsensus bandConsensus;
    TConsensus myRead;
    TConsensus newConsensus;
    int i = 0;
    for (; alignIt != alignItEnd; ++alignIt) {
        double tBegin = sysTime();
        if (i++ > 1000 || i == 1) {
            //printf("realigning %u/%u\r", unsigned(alignIt-beg), unsigned(alignItEnd-beg));
            if (i != 1)
                i = 0;
            fflush(stdout);
        }
        //// Debug code
        //for(TSize i = 0; i<length(consensus); ++i) {
        //    std::cout << consensus[i] << std::endl;
        //}
        //TAlignedReadIter debugIt = begin(contigReads, Standard() );
        //TAlignedReadIter debugItEnd = end(contigReads, Standard() );
        //for(;debugIt != debugItEnd; goNext(debugIt)) {
        //    std::cout << debugIt->beginPos << ',' << debugIt->endPos << ':';
        //    typedef typename Iterator<String<TGapAnchor> , Standard>::Type TGapIter;
        //    TGapIter gapIt = begin(debugIt->gaps, Standard());
        //    TGapIter gapItEnd = end(debugIt->gaps, Standard());
        //    for(;gapIt != gapItEnd; goNext(gapIt)) {
        //        std::cout << '(' << gapIt->seqPos << ',' << gapIt->gapPos << ')' << ',';
        //    }
        //    std::cout << std::endl;
        //}

        TSize itConsPos = 0;
        TConsIter itCons = begin(consensus, Standard());
        TConsIter itConsEnd = end(consensus, Standard());

        // Initialize the consensus of the band.
        clear(myRead);
        resize(myRead, length(fragStore.readSeqStore[alignIt->readId]), TProfileChar());
        resize(bandConsensus, 2 * bandwidth + (alignIt->endPos - alignIt->beginPos), Generous());
        TConsIter bandConsIt = begin(bandConsensus);
        TConsIter bandConsItEnd = end(bandConsensus);
        TConsIter myReadIt = begin(myRead);
        TReadPos bandOffset = 0;
        if (bandwidth < (TBandwidth) alignIt->beginPos) {
            bandOffset = alignIt->beginPos - bandwidth;
            itCons += bandOffset; itConsPos += bandOffset;
            SEQAN_ASSERT_LEQ(itCons, itConsEnd);
        }
        int leftDiag = (alignIt->beginPos - bandOffset) - bandwidth;
        int rightDiag = leftDiag + 2 * bandwidth;
        //int increaseBand = 0;
        int increaseBandLeft = 0;
        int increaseBandRight = 0;
        int removedBeginPos = 0;
        int removedEndPos = 0;
        for (TReadPos iPos = bandOffset; iPos < alignIt->beginPos && itCons != itConsEnd && bandConsIt != bandConsItEnd; ++itCons, ++bandConsIt, ++itConsPos, ++iPos)
            *bandConsIt = *itCons; // fill in positions left of readbegin
        TSize itConsPosBegin = itConsPos;  // start position of read basically, right? if(itConsPosBegin != alignIt->beginPos) std::cout <<"nicht unbedingt gleich\n";
        alignIt->beginPos = alignIt->endPos = 0; // So this read is discarded in all gap operations

        // Remove sequence from profile (and add to the consensus??)  // TODO(holtgrew): Add to consensus part right?
        typedef typename Iterator<TReadSeq, Standard>::Type TReadIter;
        TReadIter itRead = begin(fragStore.readSeqStore[alignIt->readId], Standard());
        TReadIter itReadEnd = end(fragStore.readSeqStore[alignIt->readId], Standard());
        typedef typename Iterator<String<TGapAnchor>, Standard>::Type TReadGapsIter;
        TReadGapsIter itGaps = begin(alignIt->gaps, Standard());
        TReadGapsIter itGapsEnd = end(alignIt->gaps, Standard());
        TReadPos old = 0;
        int diff = 0;
        TReadPos clippedBeginPos = 0;
        TReadPos clippedEndPos = 0;
        SEQAN_ASSERT_LT(itRead, itReadEnd);
        if ((itGaps != itGapsEnd) && (itGaps->gapPos == 0)) {
            old = itGaps->seqPos;
            clippedBeginPos = old; // gaps at beginning? or really clipped?
            //std::cout << "clippedBeginPos = " << clippedBeginPos << std::endl;
            itRead += old;
            diff -= old;
            ++itGaps;
            SEQAN_ASSERT_LT(itRead, itReadEnd);
        }
        for (; itGaps != itGapsEnd && itCons != itConsEnd; ++itGaps) {
            // limit should never be larger than read length
            TReadPos limit = itGaps->seqPos;
            SEQAN_ASSERT_LT(itGaps->seqPos, (TReadPos)length(fragStore.readSeqStore[alignIt->readId]));
            int newDiff = (itGaps->gapPos - limit);
            SEQAN_ASSERT_LT(itGaps->gapPos, (TReadPos)length(consensus));
            if (diff > newDiff) {
                clippedEndPos = diff - newDiff;
                limit -= clippedEndPos;
            }
            for (; old < limit && itCons != itConsEnd && itRead != itReadEnd && bandConsIt != bandConsItEnd; ++old, ++itRead) {
                //SEQAN_ASSERT_LT(itCons, itConsEnd);
                --(*itCons).count[ordValue(*itRead)];
                if (!empty(*itCons)) {
                    *bandConsIt = *itCons;
                    ++bandConsIt;
                    ++itConsPos;
                    removedEndPos = 0;
                } else {
                    if (itConsPosBegin != itConsPos) {
                        ++increaseBandLeft; // insertion --> increaseBandLeft, read has character here, consensus doesnt
                        ++removedEndPos;
                    } else ++removedBeginPos; // begin gaps
                    removeGap(contigReads, itConsPos);
                }
                (*myReadIt).count[0] = ordValue(*itRead);
                ++myReadIt;
                ++itCons;
                //SEQAN_ASSERT_LT(itRead, itReadEnd);
            }
            for (; diff < newDiff && itCons != itConsEnd && bandConsIt != bandConsItEnd; ++diff) {
                ++increaseBandRight; // deletion --> increaseBandRight, read has gaps here, consensus doesnt
                //SEQAN_ASSERT_LT(itCons, itConsEnd);
                --(*itCons).count[gapPos];
                if (!empty(*itCons)) {
                    *bandConsIt = *itCons;
                    ++bandConsIt;
                    ++itConsPos;
                } else removeGap(contigReads, itConsPos);  //++increaseBandRight;}
                ++itCons;
            }
        }
        if (!clippedEndPos) {
            for (; itRead!=itReadEnd && itCons != itConsEnd && bandConsIt != bandConsItEnd; ++itRead) {
                //SEQAN_ASSERT_LT(itCons, itConsEnd);
                //SEQAN_ASSERT_LT(itRead, itReadEnd);
                --(*itCons).count[ordValue(*itRead)];  //subtract the read base to get bandConsensus wo myRead
                if (!empty(*itCons)) {
                    *bandConsIt = *itCons;
                    ++bandConsIt;
                    ++itConsPos;
                    removedEndPos = 0;
                } else {  // only gaps left in this column after removing myRead
                    if (itConsPosBegin != itConsPos) {
                        ++increaseBandLeft; // insertion --> increaseBandLeft, read is longer than consensus here
                        ++removedEndPos;
                    } else ++removedBeginPos;
                    removeGap(contigReads, itConsPos);
                }
                (*myReadIt).count[0] = ordValue(*itRead);
                ++myReadIt;
                ++itCons;
            }
        }
        bool singleton = (itConsPosBegin == itConsPos);
        increaseBandLeft -= removedEndPos;
        //increaseBand = increaseBandLeft + increaseBandRight;

        // Go further up to the bandwidth
        for (TReadPos iPos = 0; ((itCons != itConsEnd) && (iPos < (TReadPos) bandwidth)) && bandConsIt != bandConsItEnd; ++itCons, ++iPos, ++bandConsIt)
            *bandConsIt = *itCons;
        resize(bandConsensus, bandConsIt - begin(bandConsensus, Standard()), Generous());
        resize(myRead, myReadIt - begin(myRead, Standard()), Generous());

        // Realign the consensus with the sequence.
        typedef StringSet<TConsensus, Dependent<> > TStringSet;
        TStringSet pairSet;
        appendValue(pairSet, bandConsensus);
        appendValue(pairSet, myRead);

        //for(TSize i = 0; i<length( pairSet[0]); ++i) {
        //    std::cout <<  pairSet[0][i] << std::endl;
        //}
        //std::cout << "_______________" << std::endl;
        //for(TSize i = 0; i<length( pairSet[1]); ++i) {
        //    std::cout <<   pairSet[1][i] << std::endl;
        //}
        //std::cout << "..............." << std::endl;

        typedef String<Fragment<> > TFragmentString;
        TFragmentString matches;
        assignProfile(consScore, bandConsensus);

        double tBegAlign = sysTime();
        leftDiag -= removedBeginPos;
        rightDiag -= removedBeginPos;
        //if(alignIt->readId == refId)
        //{
        //    std::cout << "length(Cons)=" <<  (int) length(pairSet[0]) << std::endl;
        //    std::cout << "length(Ref)=" <<  (int) length(pairSet[1]) << std::endl;
        //    std::cout << "-->leftDiag" << _max(leftDiag - increaseBandLeft, -1 * (int) length(pairSet[1])) << std::endl;
        //    std::cout << "-->rightDiag" << _min(rightDiag + increaseBandRight, (int) length(pairSet[0])) << std::endl;
        //}
        if (!singleton) {
            if (rmethod == 0)
                globalAlignment(matches, pairSet, consScore, AlignConfig<true,false,false,true>(), _max(leftDiag - increaseBandLeft, -1 * (int) length(pairSet[1])), _min(rightDiag + increaseBandRight, (int) length(pairSet[0])), NeedlemanWunsch());
            else if (rmethod == 1)
                globalAlignment(matches, pairSet, consScore, AlignConfig<true,false,false,true>(), _max(leftDiag - increaseBandLeft, -1 * (int) length(pairSet[1])), _min(rightDiag + increaseBandRight, (int) length(pairSet[0])), Gotoh());
        }
        double tEndAlign = sysTime();

        //// Debug code
        //Graph<Alignment<TStringSet, void, WithoutEdgeId> > g1(pairSet);
        //int sc1 = globalAlignment(g1, consScore, AlignConfig<true,false,false,true>(), _max(leftDiag - increaseBand, -1 * (int) length(pairSet[1])), _min(rightDiag + increaseBand, (int) length(pairSet[0])), Gotoh());
        //std::cout << sc1 << std::endl;
        //std::cout << g1 << std::endl;

        // Add the read back to the consensus and build the new consensus.
        resize(newConsensus, length(bandConsensus) + length(myRead), Generous());
        TConsIter newConsIt = begin(newConsensus, Standard());
        TConsIter bandIt = begin(bandConsensus, Standard());
        TConsIter bandItEnd = end(bandConsensus, Standard());
        typedef typename Iterator<TFragmentString, Standard>::Type TFragIter;
        TFragIter fragIt = end(matches, Standard());
        TFragIter fragItEnd = begin(matches, Standard());
        TReadPos consPos = 0;
        TReadPos readPos = 0;
        TReadPos alignPos = 0;
        clear(alignIt->gaps);
        diff = 0;
        if (clippedBeginPos) {
            appendValue(alignIt->gaps, TGapAnchor(clippedBeginPos, 0), Generous() );
            diff -= clippedBeginPos;
        }
        bool firstMatch = true;
        if (fragIt != fragItEnd) { // walk through segment matches that represent read-msa alignment
            do {
                --fragIt;
                int gapLen = fragIt->begin1 - consPos;
                if (firstMatch) gapLen = 0; // gap between two adjacent segment matches
                // equivalent to profilePos + fraglen < nextProfilePos
                while (consPos < (TReadPos)fragIt->begin1) { // cons stretch before newCons start
                    SEQAN_ASSERT_LT(bandIt, bandItEnd);
                    SEQAN_ASSERT_LT(newConsIt, end(newConsensus,Standard()));
                    if (!firstMatch) ++(*bandIt).count[gapPos]; // fill with gaps if we are between two segment matches
                    *newConsIt = *bandIt;
                    ++newConsIt;
                    ++bandIt;
                    ++consPos;
                    ++alignPos;
                }
                // equivalent to refPos + fraglen < nextRefPos
                while (readPos < (TReadPos)fragIt->begin2) { // read stretch before matching fragment starts
                    SEQAN_ASSERT_LT(readPos, (TReadPos)length(fragStore.readSeqStore[alignIt->readId]));
                    SEQAN_ASSERT_LT(newConsIt, end(newConsensus,Standard()));
                    // equivalent to profileDel
                    if (gapLen) {
                        diff += gapLen; // add gap of length gaplen to readGaps
                        appendValue(alignIt->gaps, TGapAnchor(clippedBeginPos + readPos, clippedBeginPos + readPos + diff), Generous() );
                        gapLen = 0; // do this only once
                    }
                    int numGaps = insertGap(contigReads, bandOffset + alignPos);
                    TProfileChar tmpChar;
                    ++tmpChar.count[myRead[readPos].count[0]]; // insert new column in profile
                    tmpChar.count[gapPos] += numGaps;
                    *newConsIt = tmpChar; ++newConsIt;
                    ++readPos; ++alignPos;
                }
                for (TSize i = 0; i<fragIt->len; ++i, ++bandIt, ++consPos, ++readPos, ++alignPos, ++newConsIt) {
                    SEQAN_ASSERT_LT(bandIt, bandItEnd);
                    SEQAN_ASSERT_LT(readPos, (TReadPos)length(fragStore.readSeqStore[alignIt->readId]));
                    SEQAN_ASSERT_LT(newConsIt, end(newConsensus,Standard()));
                    if (firstMatch) {
                        firstMatch = false;
                        alignIt->beginPos = bandOffset + consPos;
                    } else if (gapLen) {
                        diff += gapLen;
                        appendValue(alignIt->gaps, TGapAnchor(clippedBeginPos + readPos, clippedBeginPos + readPos + diff), Generous() );
                        gapLen = 0;
                    }
                    SEQAN_ASSERT_LT(bandIt, bandItEnd);
                    ++(*bandIt).count[myRead[readPos].count[0]];
                    *newConsIt = *bandIt;
                }
            } while (fragIt != fragItEnd);
        }

        for (; readPos < (TReadPos)length(myRead); ++readPos) {
            int numGaps = insertGap(contigReads, bandOffset + alignPos);
            TProfileChar tmpChar;
            ++tmpChar.count[myRead[readPos].count[0]];
            tmpChar.count[gapPos] += numGaps;
            SEQAN_ASSERT_LT(newConsIt, end(newConsensus,Standard()));
            *newConsIt = tmpChar; ++newConsIt;
            ++alignPos;
        }
        if (singleton) alignIt->beginPos = bandOffset;
        alignIt->endPos = alignIt->beginPos + clippedBeginPos + readPos + diff;
        if (clippedEndPos) {
            diff -= clippedEndPos;
            appendValue(alignIt->gaps, TGapAnchor(clippedBeginPos + readPos + clippedEndPos, clippedBeginPos + readPos + clippedEndPos + diff), Generous() );
        }
        for (; bandIt != bandItEnd; ++bandIt, ++newConsIt)
        {
            SEQAN_ASSERT_LT(newConsIt, end(newConsensus,Standard()));
            *newConsIt = *bandIt;
        }
        resize(newConsensus, newConsIt - begin(newConsensus, Standard()), Generous());

        replace(consensus, bandOffset, itCons - begin(consensus), newConsensus);
        double tEnd = sysTime();

        timeBeforeAlign += tBegAlign - tBegin;
        timeAlign += tEndAlign - tBegAlign;
        timeAfterAlign += tEnd -tEndAlign;
    }
}

//////////////////////////////////////////////////////////////////////////////////

/*!
 * @fn reAlign
 * @headerfile <seqan/consensus.h>
 * @brief Perform realignment using the Anson-Myers realignment.
 *
 * @deprecated Do not use this function but use the new function @link reAlignment @endlink instead.
 *
 * @signature void reAlign(fragStore, consensusScore, contigID, [realignmentMethod,] bandwidth, includeReference);
 *
 * @param[in,out] fragStore        The @link FragmentStore @endlink with the alignment to realign.
 * @param[in]     consensusScore   The @link Score @endlink to use for scoring alignments.
 * @param[in]     contigID         The integer id of the contig to realign.
 * @param[in]     bandwidth        The integer bandwidth to use for realignment.
 * @param[in]     includeReference A <tt>bool</tt> flag that indicates whether to include the reference as a pseudo-read.
 *
 * If <tt>includeReference</tt> then the reference of the given contig will be used as a pseudo-read.  In this case, the
 * reference will be replaced by the consensus.  When included as a pseudo-read, the alignment of the consensus relative
 * to the original refernence can be used to call variants.
 */

// TODO(holtgrew): realignmentMethod should not be optional or moved to the end of the list.
// TODO(holtgrew): The method should be selected with an enum instead of an int.
template<typename TSpec, typename TConfig, typename TScore, typename TId, typename TMethod, typename TBandwidth>
void
reAlign(FragmentStore<TSpec, TConfig> & fragStore,
        TScore & consScore,
        TId const contigId,
        TMethod const rmethod,
        TBandwidth const bandwidth,
        bool includeReference)
{
    typedef FragmentStore<TSpec, TConfig> TFragmentStore;
    typedef typename Size<TFragmentStore>::Type TSize;
    typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
    typedef typename TFragmentStore::TReadPos TReadPos;

    typedef typename TFragmentStore::TContigStore        TContigStore;
    typedef typename Value<TContigStore>::Type        TContig;
    typedef typename TFragmentStore::TContigPos         TContigPos;
    typedef typename TFragmentStore::TContigSeq         TContigSeq;
    typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> >    TContigGaps;

    typedef typename TFragmentStore::TReadSeq TReadSeq;
    typedef typename TFragmentStore::TReadGapAnchor TGapAnchor;
    typedef typename Value<typename TFragmentStore::TReadSeq>::Type TStoreAlphabet;
    typedef typename BaseAlphabet<TStoreAlphabet>::Type TAlphabet;
    typedef typename Value<TAlignedReadStore>::Type TAlignedElement;
    //typedef typename Value<typename TFragmentStore::TReadStore>::Type TReadElement;

    // double beginTime, endTime;

    // beginTime = sysTime();
    // TODO(holtgrew): Unnecessary, only required once.
    // Sort according to contigId
    sortAlignedReads(fragStore.alignedReadStore, SortContigId());

    typedef typename Iterator<TAlignedReadStore, Standard>::Type TAlignIter;
    TAlignIter alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
    TAlignIter alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());

    // Sort the reads according to their begin position.
    sortAlignedReads(infix(fragStore.alignedReadStore, alignIt - begin(fragStore.alignedReadStore, Standard()), alignItEnd - begin(fragStore.alignedReadStore, Standard())), SortBeginPos());
    alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
    alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
    // endTime = sysTime();
//    std::cerr << "TIME sorting " << endTime - beginTime << std::endl;

    // beginTime = sysTime();
    // Copy all reads belonging to this contig and reverse complement them if necessary.
    TAlignedReadStore contigReads;  // TODO(holtgrew): Rather contigAlignedReads?
    TReadPos maxPos = 0;
    TReadPos minPos = MaxValue<TReadPos>::VALUE;
    for (; alignIt != alignItEnd; ++alignIt) {
        if (alignIt->beginPos > alignIt->endPos) {
            reverseComplement(fragStore.readSeqStore[alignIt->readId]);
            TAlignedElement alignedEl = *alignIt;
            TReadPos tmp = alignedEl.beginPos;
            alignedEl.beginPos = alignedEl.endPos;
            alignedEl.endPos = tmp;
            if (alignedEl.beginPos < minPos)
                minPos = alignedEl.beginPos;
            if (alignedEl.endPos > maxPos)
                maxPos = alignedEl.endPos;
            appendValue(contigReads, alignedEl, Generous() );
        } else {
            if (alignIt->beginPos < minPos)
                minPos = alignIt->beginPos;
            if (alignIt->endPos > maxPos)
                maxPos = alignIt->endPos;
            appendValue(contigReads, value(alignIt), Generous() );
        }
    }

    // Append reference sequence to aligned reads for contigs if requested to do so.
    if (includeReference) {
        TId dummyReadId = length(fragStore.readSeqStore);
        TId dummyMatchId = length(fragStore.alignedReadStore);
        appendRead(fragStore, fragStore.contigStore[contigId].seq);
        appendValue(fragStore.readNameStore, fragStore.contigNameStore[contigId], Generous());
        fragStore.contigNameStore[contigId] += "Consensus_";

        TAlignedElement el;
        el.id = dummyMatchId;
        el.readId = dummyReadId;
        el.contigId = contigId;
        minPos = el.beginPos = 0;
        TContigGaps contigGaps(fragStore.contigStore[contigId].seq, fragStore.contigStore[contigId].gaps);
        maxPos = el.endPos = _max(maxPos,(TReadPos)positionSeqToGap(contigGaps,length(fragStore.contigStore[contigId].seq)-1)+1);
        maxPos = el.endPos = _max(maxPos,(TReadPos)length(fragStore.contigStore[contigId].seq));
        el.gaps = fragStore.contigStore[contigId].gaps;
        appendValue(contigReads, el, Generous());
    }
    // endTime = sysTime();
//    std::cerr << "TIME copying " << endTime - beginTime << std::endl;

    // beginTime = sysTime();
    // Create the consensus sequence
    TSize gapPos = ValueSize<TAlphabet>::VALUE;
    typedef ProfileChar<TAlphabet> TProfile;
    typedef String<TProfile> TProfileString;
    typedef typename Iterator<TProfileString, Standard>::Type TConsIter;
    TProfileString consensus;
    resize(consensus, maxPos - minPos, TProfile());

    TConsIter itCons = begin(consensus, Standard() );
    TConsIter itConsEnd = end(consensus, Standard());
    TAlignIter contigReadsIt = begin(contigReads, Standard() );
    TAlignIter contigReadsItEnd = end(contigReads, Standard() );
    for(;contigReadsIt != contigReadsItEnd; ++contigReadsIt) {
        contigReadsIt->beginPos -= minPos;
        contigReadsIt->endPos -= minPos;
        itCons = begin(consensus, Standard() );
        itCons += contigReadsIt->beginPos;

        typedef typename Iterator<TReadSeq, Standard>::Type TReadIter;
        TReadIter itRead = begin(fragStore.readSeqStore[contigReadsIt->readId], Standard() );
        TReadIter itReadEnd = end(fragStore.readSeqStore[contigReadsIt->readId], Standard() );
        typedef typename Iterator<String<typename TFragmentStore::TReadGapAnchor>, Standard>::Type TReadGapsIter;
        TReadGapsIter itGaps = begin(contigReadsIt->gaps, Standard() );
        TReadGapsIter itGapsEnd = end(contigReadsIt->gaps, Standard() );

        TReadPos old = 0;
        int diff = 0;
        bool clippedEnd = false;
        if ((itGaps != itGapsEnd) && (itGaps->gapPos == 0)) {
            old = itGaps->seqPos;
            itRead += old;
            diff -= old;
            ++itGaps;
        }
        for(;itGaps != itGapsEnd; ++itGaps) {
            TReadPos limit = itGaps->seqPos;
            int newDiff = (itGaps->gapPos - limit);
            SEQAN_ASSERT_LT(itGaps->gapPos, (int)length(consensus));
            if (diff > newDiff) {
                limit -= (diff - newDiff);
                clippedEnd = true;
            }
            for(;old < limit && itRead != itReadEnd && itCons != itConsEnd; ++old, ++itRead)
            {
                SEQAN_ASSERT_LT(itRead, itReadEnd);
                ++(value(itCons++)).count[ordValue(*itRead)];
            }
            for(;diff < newDiff; ++diff)
                ++(value(itCons++)).count[gapPos];
        }
        if (!clippedEnd) {
            for( ; itRead!=itReadEnd && itCons != itConsEnd;++itRead)
                ++(value(itCons++)).count[ordValue(*itRead)];
        }
    }
    // endTime = sysTime();
//    std::cerr << "TIME consensus " << endTime - beginTime << std::endl;

    // beginTime = sysTime();
    double tBefore = 0, tAlign = 0, tAfter = 0;
    reAlign(fragStore, contigReads, consensus, consScore, rmethod, bandwidth, includeReference, tBefore, tAlign, tAfter);
//    fprintf(stderr, "TIME before align: %f s\nTIME align: %f s\nTIME after align: %f s\n", tBefore, tAlign, tAfter);
    // endTime = sysTime();
//    std::cerr << "TIME realign " << endTime - beginTime << std::endl;
    int score = scoreConsensus(consensus);
    int oldScore = score + 1;
    while(score < oldScore) {
        //std::cout << "Score: " << score << std::endl;
        oldScore = score;
//        double beginTime = sysTime();
        double tBefore = 0, tAlign = 0, tAfter = 0;
        reAlign(fragStore, contigReads, consensus, consScore, rmethod, bandwidth, includeReference, tBefore, tAlign, tAfter);
//        fprintf(stderr, "TIME before align: %f s\nTIME align: %f s\nTIME after align: %f s\n", tBefore, tAlign, tAfter);
//        double endTime = sysTime();
//        std::cerr << "TIME realign " << endTime - beginTime << std::endl;
        score = scoreConsensus(consensus);
    }
    //std::cout << "FinalScore: " << score << std::endl;

    // beginTime = sysTime();
    // Update all the aligned reads and the new consensus
    alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
    TAlignIter contigReadIt = begin(contigReads, Standard());
    for (; alignIt != alignItEnd; ++alignIt) {
        if (alignIt->beginPos > alignIt->endPos) {
            reverseComplement(fragStore.readSeqStore[alignIt->readId]);
            alignIt->beginPos = contigReadIt->endPos;
            alignIt->endPos = contigReadIt->beginPos;
        } else {
            alignIt->beginPos = contigReadIt->beginPos;
            alignIt->endPos = contigReadIt->endPos;
        }
        // Remove empty gap anchors
        clear(alignIt->gaps);
        typedef typename Iterator<TGapAnchor, Standard>::Type TGapIter;
        TGapIter gapIt = begin(contigReadIt->gaps, Standard());
        TGapIter gapItEnd = end(contigReadIt->gaps, Standard());
        int diff = 0;
        for(;gapIt != gapItEnd; ++gapIt) {
            if ((int) gapIt->gapPos - (int) gapIt->seqPos != diff) {
                diff = (int) gapIt->gapPos - (int) gapIt->seqPos;
                appendValue(alignIt->gaps, *gapIt, Generous() );
            }
        }
        ++contigReadIt;
    }
    typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
    TContigElement& contigEl = fragStore.contigStore[contigId];
    typedef typename Iterator<TProfileString, Standard>::Type TConsIter;
    TConsIter itConsensus = begin(consensus, Standard());
    TConsIter itConsensusEnd = end(consensus, Standard());
    char gapChar = gapValue<char>();
    TSize gapLen = 0;
    TContigPos contigPos = 0;
    int diff = 0;
    clear(contigEl.seq);
    clear(contigEl.gaps);
    for (; itConsensus != itConsensusEnd; ++itConsensus, ++contigPos) {
        if ((char) *itConsensus == gapChar) ++gapLen;
        else {
            if (gapLen) {
                diff += (int) gapLen;
                appendValue(contigEl.gaps, TGapAnchor(contigPos - diff, contigPos), Generous() );
                gapLen = 0;
            }
            // TODO(weese): Here we convert from ProfileChar<Dna5>->Dna5->Dna5Q
            // instead diverting through Dna5 we could think of directly converting
            // a profile to a quality value, e.g. like the base caller Phred does.
            // Therefore the conversion ProfileChar<Dna5> <-> Dna5Q needs to be
            // defined.
            appendValue(contigEl.seq, (TAlphabet)value(itConsensus), Generous() );
        }

    }
    if (includeReference)
        appendValue(fragStore.alignedReadStore, contigReads[length(contigReads) - 1], Generous() );
    // endTime = sysTime();
//    std::cerr << "TIME finalizing " << endTime - beginTime << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////

// Forwards to the overload that accepts the alignment method.
template<typename TSpec, typename TConfig, typename TScore, typename TId, typename TBandwidth>
inline void
reAlign(FragmentStore<TSpec, TConfig>& fragStore,
        TScore& consScore,
        TId const contigId,
        TBandwidth const bandwidth,
        bool includeReference)
{
    reAlign(fragStore, consScore, contigId, 0, bandwidth, includeReference);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CONSENSUS_CONSENSUS_REALIGN_H_
