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

#ifndef SEQAN_HEADER_CONSENSUS_LIBRARY_H
#define SEQAN_HEADER_CONSENSUS_LIBRARY_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

template<typename TSize, typename TCargo1, typename TCargo2>
inline void
_getAlignmentStatistics(Nothing&,
                         TSize,
                         TSize,
                         TSize,
                         TCargo1,
                         TCargo2)
{
    SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSpec, typename TSize, typename TCargo1, typename TCargo2>
inline void
_getAlignmentStatistics(String<TValue, TSpec>& dist,
                         TSize i,
                         TSize j,
                         TSize nseq,
                         TCargo1,
                         TCargo2 quality)
{
    SEQAN_CHECKPOINT
    dist[i*nseq + j] = (TValue) (100 - quality);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TSize, typename TCargo1, typename TCargo2>
inline void
_getAlignmentStatistics(Graph<Undirected<TCargo, TSpec> >& dist,
                         TSize i,
                         TSize j,
                         TSize,
                         TCargo1,
                         TCargo2 quality)
{
    SEQAN_CHECKPOINT
    addEdge(dist, i, j, (TCargo) (100 - quality));
}

//////////////////////////////////////////////////////////////////////////////
// Layout-based pair selection
//////////////////////////////////////////////////////////////////////////////

template<typename TSize>
struct LessPair_ :
    public std::unary_function<Pair<TSize, TSize>, bool>
{
    inline bool
    operator() (Pair<TSize, TSize> const& a1, Pair<TSize, TSize> const& a2) const {
        if (a1.i1 == a2.i1) return (a1.i2 < a2.i2);
        else return (a1.i1 < a2.i1);
    }
};

template<typename TSize>
struct _LessTripel :
    public std::unary_function<Pair<TSize, Triple<TSize, TSize, TSize> >, bool>
{
    inline bool
    operator() (Pair<TSize, Triple<TSize, TSize, TSize> > const& a1, Pair<TSize, Triple<TSize, TSize, TSize> > const& a2) {
        return (a1.i1 < a2.i1);
    }
};


//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TBegEndPos, typename TSize, typename TPairList, typename TPos, typename TSpec2>
inline void
selectPairsAssembly(StringSet<TString, TSpec> const & str,
                   TBegEndPos const & begEndPos,
                   TSize bandwidth,
                   TPairList & pList,
                   String<Pair<TPos, TPos>, TSpec2> & dList)
{
    typedef String<Pair<TPos, TPos>, TSpec2>  TDistanceList;
    //typedef StringSet<TString, TSpec> TStringSet;
    typedef Pair<TPos, TPos> TDiagPair;
    typedef typename Value<TPairList>::Type TPair;
    //typedef typename Iterator<TPairList, Standard>::Type TPairIter;
    typedef typename Iterator<TBegEndPos const, Standard>::Type TBegEndIter;

    // Initialization
    TSize nseq = length(str);

    // Workaround for strange celera behaviour (just for contained reads)
#ifdef CELERA_OFFSET
    TSize contained_offset=200;
#else
    TSize contained_offset=0;
#endif

    // Sort the reads by their first index position
    TBegEndIter begEndIt = begin(begEndPos, Standard());
    TBegEndIter begEndItEnd = end(begEndPos, Standard());
    typedef Triple<TSize, TSize, TSize> TInfo;
    typedef String<Pair<TSize, TInfo> > TPosIndexList;
    typedef typename Iterator<TPosIndexList, Standard>::Type TPosIter;
    TPosIndexList posIndex;
    resize(posIndex, length(begEndPos));
    TPosIter posIndexIt = begin(posIndex, Standard());
    TPosIter posIndexItEnd = end(posIndex, Standard());
    for(TSize index = 0;begEndIt != begEndItEnd; ++begEndIt, ++posIndexIt, ++index)
        *posIndexIt = ((*begEndIt).i1 < (*begEndIt).i2) ? Pair<TSize, TInfo>((*begEndIt).i1,TInfo(index, (*begEndIt).i1, (*begEndIt).i2)) : Pair<TSize, TInfo>((*begEndIt).i2,TInfo(index, (*begEndIt).i1, (*begEndIt).i2));
    std::sort(begin(posIndex, Standard() ), end(posIndex, Standard() ), _LessTripel<TSize>() );

    // The expected overlap by a pair of reads (represented by its index)
    typedef String<Pair<TSize, TSize> > TOverlapIndexList;
    TOverlapIndexList ovlIndex;
    TSize pairLen = 0;  // Pair Counter
    TPos const initialRadius = (bandwidth + 1) / 2;
    TPos const lengthDivider = 5;    // Overlap / 2^lengthDivider is added to the radius

    // Find all overlapping reads
    TDistanceList preDList;
    TPairList prePList;
    reserve(preDList, nseq * 40);
    reserve(prePList, nseq * 40);
    reserve(ovlIndex, nseq * 40);
    posIndexIt = begin(posIndex, Standard());
    posIndexItEnd = end(posIndex, Standard());
    for(;posIndexIt != posIndexItEnd; ++posIndexIt) {
        TSize index1 = ((*posIndexIt).i2).i1;
        TPos posIi1 = ((*posIndexIt).i2).i2;
        TPos posIi2 = ((*posIndexIt).i2).i3;
        bool forwardI = (posIi1 < posIi2) ? true : false;
        TSize lenI = (posIi1 < posIi2) ? posIi2 - posIi1 : posIi1 - posIi2;
        TPosIter posIndexIt2 = posIndexIt;
        ++posIndexIt2;
        for(;posIndexIt2 != posIndexItEnd; ++posIndexIt2) {
            if ((*posIndexIt).i1 + lenI <= (*posIndexIt2).i1) break;
            TSize index2 = ((*posIndexIt2).i2).i1;
            TPos posJi1 = ((*posIndexIt2).i2).i2;
            TPos posJi2 = ((*posIndexIt2).i2).i3;

            // Diagonal boundaries of the band
            // Initialization values are used if one read is contained in the other
            TSize lenJ = (posJi1 < posJi2) ? posJi2 - posJi1 : posJi1 - posJi2;
            bool forwardJ = (posJi1 < posJi2) ? true : false;
            TPos diagLow = -1 * (TPos) lenJ;
            TPos diagHigh = (TPos) lenI;
            TPos radius = initialRadius;        // Increased by overlap length

            // Read orientations
            if (forwardI) {
                // 1) Forward - Forward
                if (forwardJ) {
                    if ((posJi2 < posIi2) && (posJi1 < posIi1)) {
                        TPos offset = (posIi1 - posJi1);
                        radius += (posJi2 - posJi1 - offset) >> lengthDivider;
                        if (-1 * offset - radius > diagLow) diagLow = -1 * offset - radius;
                        if (-1 * offset + radius < diagHigh) diagHigh = -1 * offset + radius;
                    } else if ((posJi1 > posIi1) && (posJi2 > posIi2)) {
                        TPos offset = (posJi1 - posIi1);
                        radius += (posIi2 - posIi1 - offset) >> lengthDivider;
                        if (offset + radius < diagHigh) diagHigh = offset + radius;
                        if (offset - radius > diagLow) diagLow = offset - radius;
                    } else {
                        radius += contained_offset;
                        if (posIi1 < posJi1) {
                            TPos offset = (posJi1 - posIi1);
                            radius += (posJi2 - posJi1) >> lengthDivider;
                            if (offset + radius < diagHigh) diagHigh = offset + radius;
                            if (offset - radius > diagLow) diagLow = offset - radius;
                        } else {
                            TPos offset = (posIi1 - posJi1);
                            radius += (posIi2 - posIi1) >> lengthDivider;
                            if (-1 * offset + radius < diagHigh) diagHigh = -1 * offset + radius;
                            if (-1 * offset - radius > diagLow) diagLow = -1 * offset - radius;
                        }
                    }
                } else { // 2) Forward - Reverse
                    if ((posJi1 < posIi2) && (posJi2 < posIi1)) {
                        TPos offset = (posIi1 - posJi2);
                        radius += (posJi1 - posJi2 - offset) >> lengthDivider;
                        if (-1 * offset + radius < diagHigh) diagHigh = -1 * offset + radius;
                        if (-1 * offset - radius > diagLow) diagLow = -1 * offset - radius;
                    } else if ((posJi2 > posIi1) && (posJi1 > posIi2)) {
                        TPos offset = (posJi2 - posIi1);
                        radius += (posIi2 - posIi1 - offset) >> lengthDivider;
                        if (offset + radius < diagHigh) diagHigh = offset + radius;
                        if (offset - radius > diagLow) diagLow = offset - radius;
                    } else {
                        radius += contained_offset;
                        if (posIi1 < posJi2) {
                            TPos offset = (posJi2 - posIi1);
                            radius += (posJi1 - posJi2) >> lengthDivider;
                            if (offset + radius < diagHigh) diagHigh = offset + radius;
                            if (offset - radius > diagLow) diagLow = offset - radius;
                        } else {
                            TPos offset = (posIi1 - posJi2);
                            radius += (posIi2 - posIi1) >> lengthDivider;
                            if (-1 * offset + radius < diagHigh) diagHigh = -1 * offset + radius;
                            if (-1 * offset - radius > diagLow) diagLow = -1 * offset - radius;
                        }
                    }
                }
            } else {
                // 3) Reverse - Forward
                if (forwardJ) {
                    if ((posIi1 > posJi2) && (posIi2 > posJi1)) {
                        TPos offset = (posIi2 - posJi1);
                        radius += (posJi2 - posJi1 - offset) >> lengthDivider;
                        if (-1 * offset - radius > diagLow) diagLow = -1 * offset - radius;
                        if (-1 * offset + radius < diagHigh) diagHigh = -1 * offset + radius;
                    } else if ((posJi1 > posIi2) && (posJi2 > posIi1)) {
                        TPos offset = (posJi1 - posIi2);
                        radius += (posIi1 - posIi2 - offset) >> lengthDivider;
                        if (offset + radius < diagHigh) diagHigh = offset + radius;
                        if (offset - radius > diagLow) diagLow = offset - radius;
                    } else {
                        radius += contained_offset;
                        if (posIi2 < posJi1) {
                            TPos offset = (posJi1 - posIi2);
                            radius += (posJi2 - posJi1) >> lengthDivider;
                            if (offset + radius < diagHigh) diagHigh = offset + radius;
                            if (offset - radius > diagLow) diagLow = offset - radius;
                        } else {
                            TPos offset = (posIi2 - posJi1);
                            radius += (posIi1 - posIi2) >> lengthDivider;
                            if (-1 * offset + radius < diagHigh) diagHigh = -1 * offset + radius;
                            if (-1 * offset - radius > diagLow) diagLow = -1 * offset - radius;
                        }
                    }
                } else { // 4) Reverse - Reverse
                    if ((posJi1 < posIi1) && (posJi2 < posIi2)) {
                        TPos offset = (posIi2 - posJi2);
                        radius += (posJi1 - posJi2 - offset) >> lengthDivider;
                        if (-1 * offset + radius < diagHigh) diagHigh = -1 * offset + radius;
                        if (-1 * offset - radius > diagLow) diagLow = -1 * offset - radius;
                    } else if ((posJi2 > posIi2) && (posJi1 > posIi1)) {
                        TPos offset = (posJi2 - posIi2);
                        radius += (posIi1 - posIi2 - offset) >> lengthDivider;
                        if (offset + radius < diagHigh) diagHigh = offset + radius;
                        if (offset - radius > diagLow) diagLow = offset - radius;
                    } else {
                        radius += contained_offset;
                        if (posIi2 < posJi2) {
                            TPos offset = (posJi2 - posIi2);
                            radius += (posJi1 - posJi2) >> lengthDivider;
                            if (offset + radius < diagHigh) diagHigh = offset + radius;
                            if (offset - radius > diagLow) diagLow = offset - radius;
                        } else {
                            TPos offset = (posIi2 - posJi2);
                            radius += (posIi1 - posIi2) >> lengthDivider;
                            if (-1 * offset + radius < diagHigh) diagHigh = -1 * offset + radius;
                            if (-1 * offset - radius > diagLow) diagLow = -1 * offset - radius;
                        }
                    }
                }
            }

            // Append this pair of reads
            if (index1 < index2) {
                appendValue(prePList, TPair(positionToId(str, index1), positionToId(str, index2)), Generous());
                appendValue(preDList, TDiagPair(diagLow, diagHigh), Generous());
            } else {
                appendValue(prePList, TPair(positionToId(str, index2), positionToId(str, index1)), Generous());
                appendValue(preDList, TDiagPair(-1 * diagHigh, -1 * diagLow), Generous());
            }
            // Estimate the overlap quality
            TPos avgDiag = (diagLow + diagHigh) / 2;
            if (avgDiag < 0) avgDiag *= -1;
            appendValue(ovlIndex, Pair<TSize, TSize>((TSize) (avgDiag), pairLen), Generous());
            ++pairLen;
        }
    }

    // Sort the pairs, better expected overlaps come first
    std::sort(begin(ovlIndex, Standard() ), end(ovlIndex, Standard() ), LessPair_<TSize>() );
    typedef typename Iterator<TOverlapIndexList, Standard>::Type TOVLIter;
    TOVLIter itOvl = begin(ovlIndex, Standard());
    TOVLIter itOvlEnd = end(ovlIndex, Standard());
    reserve(dList, pairLen);
    reserve(pList, pairLen);
    for(;itOvl != itOvlEnd; ++itOvl) {
        TSize count = (*itOvl).i2;
        appendValue(dList, preDList[count], Generous());
        appendValue(pList, prePList[count], Generous());
    }
}



//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TBegEndPos, typename TSize, typename TPairList, typename TPos, typename TSpec2>
inline void
selectPairsAllAgainstAll(StringSet<TString, TSpec> const & str,
                         TBegEndPos const & begEndPos,
                         TSize lookAround,
                         TPairList& pList,
                         String<Pair<TPos, TPos>, TSpec2> & dList)
{
    //typedef String<Pair<TPos, TPos>, TSpec2>  TDistanceList;
    //typedef StringSet<TString, TSpec> TStringSet;
    typedef Pair<TPos, TPos> TDiagPair;
    typedef typename Value<TPairList>::Type TPair;
    typedef typename Iterator<TBegEndPos const, Standard>::Type TBegEndIter;

    TBegEndIter beIt = begin(begEndPos, Standard());
    TBegEndIter beItEnd = end(begEndPos, Standard());
    TSize index1 = 0;
    for(;beIt != beItEnd; ++beIt, ++index1) {
        TPos beg2 = (beIt->i2 < beIt->i1) ? beIt->i2 : beIt->i1;
        TPos diagHigh = (beIt->i2 < beIt->i1) ? beIt->i1 - beIt->i2 : beIt->i2 - beIt->i1;
        TBegEndIter beIt2 = beIt;
        TSize index2 = index1 + 1;
        for(;++beIt2 != beItEnd; ++index2) {
            TPos beg1 = (beIt2->i2 < beIt2->i1) ? beIt2->i2 : beIt2->i1;
            TPos diagLow = (beIt2->i2 < beIt2->i1) ? beIt2->i2 - beIt2->i1 : beIt2->i1 - beIt2->i2;
            TPos diff = (beg1 > beg2) ? beg1 - beg2 : beg2 - beg1;
            if (diff < (TPos) lookAround) {
                appendValue(pList, TPair(positionToId(str, index1), positionToId(str, index2)), Generous());
                appendValue(dList, TDiagPair(diagLow, diagHigh), Generous());
            }
        }
    }
}


//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TId, typename TDiagList, typename TBegEndPos, typename TScore, typename TSize, typename TSegmentMatches, typename TScoreValues, typename TDistance>
inline void
appendSegmentMatches(StringSet<TString, TSpec> const & str,
                     String<Pair<TId, TId> > const & pList,
                     TDiagList const & dList,
                     TBegEndPos const & begEndPos,
                     TScore const & score_type,
                     TSize thresholdMatchlength,
                     TSize thresholdQuality,
                     TSize maxOvl,
                     TSegmentMatches & matches,
                     TScoreValues & scores,
                     TDistance & dist,
                     OverlapLibrary)
{
    typedef StringSet<TString, Dependent<> > TStringSet;
    typedef String<Pair<TId, TId> > const TPairList;
    typedef typename Value<TScoreValues>::Type TScoreValue;
    typedef typename Iterator<TPairList, Standard>::Type TPairIter;
    typedef typename Iterator<TDiagList const, Standard>::Type TDiagIter;

    // Initialization
    TSize nseq = length(str);
    _resizeWithRespectToDistance(dist, nseq);

    // "Front" and "Back"-overlap counter for each read
    String<TSize> frontOvl;
    String<TSize> backOvl;
    resize(frontOvl, nseq, 0);
    resize(backOvl, nseq, 0);

    // Pairwise alignments
    String<bool> aligned;
    resize(aligned, length(pList), true);
    typedef Iterator<String<bool>, Standard>::Type TBoolIter;
    TBoolIter itAligned = begin(aligned, Standard());
    TPairIter itPair = begin(pList, Standard());
    TDiagIter itDiag = begin(dList, Standard());
    TPairIter itPairEnd = end(pList, Standard());
    TSize dropCount = 0;
    for(;itPair != itPairEnd; ++itPair, ++itDiag, ++itAligned) {
        TId id1 = itPair->i1;
        TId id2 = itPair->i2;
        TSize seq1 = idToPosition(str, id1);
        TSize seq2 = idToPosition(str, id2);
        if ((frontOvl[seq1] > maxOvl) && (backOvl[seq1] > maxOvl) &&
            (frontOvl[seq2] > maxOvl) && (backOvl[seq2] > maxOvl)) {
                ++dropCount;
                continue;
        }


        // Make a pairwise string-set
        TStringSet pairSet;
        assignValueById(pairSet, const_cast<StringSet<TString, TSpec>&>(str), id1);
        assignValueById(pairSet, const_cast<StringSet<TString, TSpec>&>(str), id2);


        // Overlap alignment
        TSize from = length(matches);
        TScoreValue myScore = globalAlignment(matches, pairSet, score_type, AlignConfig<true,true,true,true>(), itDiag->i1, itDiag->i2, Gotoh() );
        TSize to = length(matches);

        // Determine a sequence weight
        TSize matchLen = 0;
        TSize overlapLen = 0;
        TSize alignLen = 0;
        if (from != to) getAlignmentStatistics(matches, pairSet, from, to, matchLen, overlapLen, alignLen);


        // Get only the good overlap alignments
        if ((overlapLen) && ((matchLen * 100) / overlapLen >= thresholdQuality) && (matchLen >= thresholdMatchlength)) {

            //// Debug Code
            //Graph<Alignment<TStringSet, TSize> > tmp(pairSet);
            //globalAlignment(tmp, pairSet, score_type, AlignConfig<true,true,true,true>(), (value(itDiag)).i1, (value(itDiag)).i2, Gotoh() );
            ////globalAlignment(tmp, pairSet, score_type, Gotoh() );
            //std::cout << "Match length: " << matchLen << std::endl;
            //std::cout << "Overlap length: " << overlapLen << std::endl;
            //std::cout << "Align length: " << alignLen << std::endl;
            //std::cout << "Quality: " << quality << std::endl;
            //std::cout << tmp << std::endl;

            // Create a corresponding edge
            if (seq1<seq2) _getAlignmentStatistics(dist, seq1, seq2, nseq, matchLen, (matchLen * 100) / overlapLen);
            else _getAlignmentStatistics(dist, seq2, seq1, nseq, matchLen, (matchLen * 100) / overlapLen);

            // Record the scores
            resize(scores, to);
            typedef typename Iterator<TScoreValues, Standard>::Type TScoreIter;
            TScoreIter itScore = begin(scores, Standard());
            TScoreIter itScoreEnd = end(scores, Standard());
            itScore += from;
            for(;itScore != itScoreEnd; ++itScore)
                *itScore = myScore;

            // Update the overlap counter
            TSize lenLast = matches[from].len;
            if (matches[to - 1].begin1 == 0) ++frontOvl[seq1];
            if (matches[to - 1].begin2 == 0) ++frontOvl[seq2];
            if (matches[from].begin1 + lenLast == length(pairSet[0])) ++backOvl[seq1];
            if (matches[from].begin2 + lenLast == length(pairSet[1])) ++backOvl[seq2];
        } else {
            resize(matches, from);
            *itAligned = false;
        }
    }
    //std::cout << "Filtration ration: " << (double) dropCount / (double) length(pList) << std::endl;

    // Find sequences that have no overlap in the front or back
    String<TSize> noFront;
    String<TSize> noBack;
    for(TSize seqI = 0; seqI < nseq; ++seqI) {
        if (frontOvl[seqI] == 0) appendValue(noFront, seqI, Generous());
        else if (backOvl[seqI] == 0) appendValue(noBack, seqI, Generous());
    }
    // Drop the first and the last sequence
    typedef typename Iterator<TBegEndPos const, Standard>::Type TBegEndIter;
    TBegEndIter begEndIt = begin(begEndPos, Standard());
    TBegEndIter begEndItEnd = end(begEndPos, Standard());
    TSize minVal = maxValue<TSize>();
    TSize maxVal = 0;
    for(;begEndIt != begEndItEnd; ++begEndIt) {
        TSize pos1 = begEndIt->i1;
        TSize pos2 = begEndIt->i2;
        if (pos1 > pos2) { TSize tmp = pos1; pos1 = pos2; pos2 = tmp;}
        if (pos1 < minVal) minVal = pos1;
        if (pos2 > maxVal) maxVal = pos2;
    }
    // Insert all remaining sequences into a string
    String<TSize> unalignedReads;
    for(TSize i = 0; i < (TSize) length(noFront); ++i) {
        TSize p1 = begEndPos[noFront[i]].i1;
        TSize p2 = begEndPos[noFront[i]].i2;
        if (p1 > p2) {TSize tmp = p1; p1 = p2; p2 = tmp; }
        if (p1 != minVal) appendValue(unalignedReads, noFront[i], Generous());
    }
    for(TSize i = 0; i < (TSize) length(noBack); ++i) {
        TSize p1 = begEndPos[noBack[i]].i1;
        TSize p2 = begEndPos[noBack[i]].i2;
        if (p1 > p2) {TSize tmp = p1; p1 = p2; p2 = tmp; }
        if (p2 != maxVal) appendValue(unalignedReads, noBack[i], Generous());
    }
    TSize countUnalignedReads = length(unalignedReads);
    //std::cout << "Unaligned reads: " << countUnalignedReads << std::endl;
    if (countUnalignedReads > 0) {
        // Sort unaligned reads
        std::sort(begin(unalignedReads, Standard()), end(unalignedReads, Standard()));

        // Realign all unaligned sequences
        itPair = begin(pList, Standard());
        itDiag = begin(dList, Standard());
        itAligned = begin(aligned, Standard());
        for(;itPair != itPairEnd; ++itPair, ++itDiag, ++itAligned) {
            if (*itAligned == true) continue;
            TId id1 = itPair->i1;
            TId id2 = itPair->i2;
            TSize seq1 = idToPosition(str, id1);
            TSize seq2 = idToPosition(str, id2);
            if ((!std::binary_search(begin(unalignedReads, Standard()), end(unalignedReads, Standard()), seq1)) &&
                (!std::binary_search(begin(unalignedReads, Standard()), end(unalignedReads, Standard()), seq2)))
                continue;
            if ((frontOvl[seq1] > maxOvl) && (backOvl[seq1] > maxOvl) &&
                (frontOvl[seq2] > maxOvl) && (backOvl[seq2] > maxOvl))
                continue;

            // Make a pairwise string-set
            TStringSet pairSet;
            assignValueById(pairSet, const_cast<StringSet<TString, TSpec>&>(str), id1);
            assignValueById(pairSet, const_cast<StringSet<TString, TSpec>&>(str), id2);

            // Overlap alignment
            TSize from = length(matches);
            TScoreValue myScore = globalAlignment(matches, pairSet, score_type, AlignConfig<true,true,true,true>(), Gotoh() );
            TSize to = length(matches);

            // Determine a sequence weight
            TSize matchLen = 0;
            TSize overlapLen = 0;
            TSize alignLen = 0;
            getAlignmentStatistics(matches, pairSet, from, to, matchLen, overlapLen, alignLen);

            if (((matchLen * 100) / overlapLen >= 80) && (matchLen >= 5)) {
                // Create a corresponding edge
                if (seq1<seq2) _getAlignmentStatistics(dist, seq1, seq2, nseq, matchLen, (matchLen * 100) / overlapLen);
                else _getAlignmentStatistics(dist, seq2, seq1, nseq, matchLen, (matchLen * 100) / overlapLen);

                // Record the scores
                resize(scores, to);
                typedef typename Iterator<TScoreValues, Standard>::Type TScoreIter;
                TScoreIter itScore = begin(scores, Standard());
                TScoreIter itScoreEnd = end(scores, Standard());
                itScore+=from;
                for(;itScore != itScoreEnd; ++itScore) *itScore = myScore;

                // Update the overlap counter
                TSize lenLast = matches[from].len;
                if (matches[to - 1].begin1 == 0) ++frontOvl[seq1];
                if (matches[to - 1].begin2 == 0) ++frontOvl[seq2];
                if (matches[from].begin1 + lenLast == length(pairSet[0])) ++backOvl[seq1];
                if (matches[from].begin2 + lenLast == length(pairSet[1])) ++backOvl[seq2];
            } else resize(matches, from);
        }
    }
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
