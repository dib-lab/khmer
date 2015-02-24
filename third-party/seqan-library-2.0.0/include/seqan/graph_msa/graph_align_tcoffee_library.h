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

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_LIBRARY_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_LIBRARY_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Alignment graph generation
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Tags
//////////////////////////////////////////////////////////////////////////////

/*!
 * @defgroup SegmentMatchGenerationTags Segment Match Generation Tags
 * @brief Tags specifying how to generate segment matches.
 *
 *
 * @tag SegmentMatchGenerationTags#GlobalPairwiseLibrary
 * @headerfile <seqan/graph_msa.h>
 * @brief Segment matches from pairwise global alignment.
 *
 * @signature typedef Tag<GlobalPairwiseLibrary_> const GlobalPairwiseLibrary;
 *
 *
 * @tag SegmentMatchGenerationTags#LocalPairwiseLibrary
 * @headerfile <seqan/graph_msa.h>
 * @brief Segment matches from pairwise local alignment.
 *
 *
 * @tag SegmentMatchGenerationTags#KmerLibrary
 * @headerfile <seqan/graph_msa.h>
 * @brief Segment matches from pairwise k-mer library.
 *
 * @signature typedef Tag<KmerLibrary_> const KmerLibrary;
 *
 *
 * @tag SegmentMatchGenerationTags#LcsLibrary
 * @headerfile <seqan/graph_msa.h>
 * @brief Segment matches from pairwise longest common subsequence comparison.
 *
 * @signature typedef Tag<LcsLibrary_> const LcsLibrary;
 */


struct GlobalPairwiseLibrary_;
typedef Tag<GlobalPairwiseLibrary_> const GlobalPairwiseLibrary;


struct LocalPairwiseLibrary_;
typedef Tag<LocalPairwiseLibrary_> const LocalPairwiseLibrary;

struct KmerLibrary_;
typedef Tag<KmerLibrary_> const KmerLibrary;


struct LcsLibrary_;
typedef Tag<LcsLibrary_> const LcsLibrary;

// ---------------------------------------------------------------------------
// BEGIN OF TO BE REMOVED LEGACY CODE
// ---------------------------------------------------------------------------

// TODO(holtgrew): Reproduction from old graph_align module. This should go away once we can compute and enumerate local matches into Fragment strings.

template<typename TSpec, typename TSize>
inline bool
_isClumping(String<bool, TSpec> const& forbidden,
            TSize row,
            TSize col,
            TSize len2)
{
    return forbidden[(col-1) * len2 + (row-1)];
}

template<typename TSize>
inline bool
_isClumping(Nothing&,
            TSize,
            TSize,
            TSize)
{
    return false;
}

struct SmithWatermanClump_;
typedef Tag<SmithWatermanClump_> SmithWatermanClump;

template <typename TAlign, typename TStringSet, typename TTrace, typename TVal, typename TIndexPair, typename TForbidden>
inline void
_alignSmithWatermanTrace(TAlign& align,
                            TStringSet const& str,
                            TTrace const& trace,
                            TVal const initialDir,
                            TIndexPair const& indexPair,
                            TForbidden& forbidden)
{
    SEQAN_CHECKPOINT
    typedef typename Size<TTrace>::Type TSize;
    typedef typename Value<TTrace>::Type TTraceValue;
    typedef typename Id<TStringSet>::Type TId;

    // TraceBack values for Gotoh
    TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2; TTraceValue Stop = 3;

    TId id1 = positionToId(const_cast<TStringSet&>(str), 0);
    TId id2 = positionToId(const_cast<TStringSet&>(str), 1);
    TSize len1 = indexPair[1];
    TSize len2 = indexPair[0];
    if ((indexPair[0] == 0) || (indexPair[1] == 0)) return;
    TSize numCols = length(str[0]);
    TSize numRowsOrig = length(str[1]);
    if (len1 < numCols) _alignTracePrint(align, str[0], str[1], id1, len1, id2, len2, numCols - len1, Horizontal);
    if (len2 < numRowsOrig) _alignTracePrint(align, str[0], str[1], id1, len1, id2, len2, numRowsOrig - len2, Vertical);
    TSize numRows = (numRowsOrig >> 1) + (numRowsOrig & 1);



    // Initialize everything
    TTraceValue nextTraceValue = (len2 & 1) ? trace[(len1 - 1)*numRows + ((len2 - 1) >> 1)] >> 4 : trace[(len1 - 1)*numRows + ((len2 - 1) >> 1)];
    TTraceValue tv = Diagonal;
    if (initialDir == Diagonal) tv = (nextTraceValue & 3);
    else if (initialDir == Horizontal) {
        if ((nextTraceValue >> 2) & 1) _alignTracePrint(align, str[0], str[1], id1, --len1, id2, len2, (TSize) 1, Horizontal);
        else tv = Horizontal;
    } else if (initialDir == Vertical) {
        if ((nextTraceValue >> 3) & 1) _alignTracePrint(align, str[0], str[1], id1, len1, id2, --len2, (TSize) 1, Vertical);
        else tv = Vertical;
    }
    TSize segLen = 0;
    TTraceValue tvOld = tv;

    // Now follow the trace
    do {
        nextTraceValue = (len2 & 1) ? trace[(len1 - 1)*numRows + ((len2 - 1) >> 1)] >> 4 : trace[(len1 - 1)*numRows + ((len2 - 1) >> 1)];
        if ((nextTraceValue & 3) == Stop) break;
        _setForbiddenCell(forbidden, len1, len2, numRowsOrig);
        if (tv == Diagonal) tv = (nextTraceValue & 3);
        else if (tv == Horizontal) {
            if ((nextTraceValue >> 2) & 1) tv = Diagonal;
            else tv =  Horizontal;
        } else if (tv == Vertical) {
            if ((nextTraceValue >> 3) & 1) tv =  Diagonal;
            else tv =  Vertical;
        }
        if (tv == Diagonal) {
            if (tv != tvOld) {
                if (tvOld == Vertical) --len2;
                else --len1;
                _alignTracePrint(align, str[0], str[1], id1, len1, id2, len2, ++segLen, tvOld);
                tvOld = tv; segLen = 0;
            } else {
                ++segLen;
                --len1; --len2;
            }
        } else if (tv == Horizontal) {
            if (tv != tvOld) {
                _alignTracePrint(align, str[0], str[1], id1, len1, id2, len2, segLen, tvOld);
                if ((nextTraceValue >> 2) & 1) {
                    _alignTracePrint(align, str[0], str[1], id1, --len1, id2, len2, (TSize) 1, Horizontal);
                    tv = Diagonal; segLen = 0;
                } else {
                    tvOld = tv; segLen = 1;
                    --len1;
                }
            } else {
                ++segLen;
                --len1;
            }
        } else if (tv == Vertical) {
            if (tv != tvOld) {
                _alignTracePrint(align, str[0], str[1], id1, len1, id2, len2, segLen, tvOld);
                if ((nextTraceValue >> 3) & 1) {
                    _alignTracePrint(align, str[0], str[1], id1, len1, id2, --len2, (TSize) 1, Vertical);
                    tv = Diagonal; segLen = 0;
                } else {
                    tvOld = tv; segLen = 1;
                    --len2;
                }
            } else {
                ++segLen;
                --len2;
            }
        }
    } while ((len1 != 0) && (len2 !=0));
    // Process left-overs
    if (segLen) _alignTracePrint(align, str[0], str[1], id1, len1, id2, len2, segLen, tvOld);

    // Handle the remaining sequence
    if (len1 != 0) _alignTracePrint(align, str[0], str[1], (TId) id1, (TSize) 0, (TId) 0, (TSize) 0, (TSize) len1, Horizontal);
    if (len2 != 0) _alignTracePrint(align, str[0], str[1], (TId) 0, (TSize) 0, (TId) id2, (TSize) 0, (TSize) len2, Vertical);
}

template <typename TTrace, typename TStringSet, typename TScore, typename TIndexPair, typename TForbidden>
inline typename Value<TScore>::Type
_alignSmithWaterman(TTrace& trace,
                      TStringSet const& str,
                      TScore const & sc,
                      typename Value<TTrace>::Type& initialDir,
                      TIndexPair& indexPair,
                      TForbidden& forbidden)
{
    SEQAN_CHECKPOINT
    typedef typename Size<TTrace>::Type TSize;
    typedef typename Value<TTrace>::Type TTraceValue;

    // TraceBack values for Smith Waterman
    TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2; TTraceValue Stop = 3;

    // The DP Matrix for diagonal walks
    typedef typename Value<TScore>::Type TScoreValue;
    typedef String<TScoreValue> TColumn;
    TColumn mat;
    // The DP Matrix for gaps from the left
    TColumn horizontal;
    // The DP Matrix for gaps from the top
    TScoreValue vert = 0;

    typedef typename Iterator<TColumn, Standard>::Type TMatIter;

    // Initialization
    typedef typename Value<TStringSet>::Type TString;
    TString const& str1 = str[0];
    TString const& str2 = str[1];
    TSize len1 = length(str1);
    TSize len2 = length(str2);
    resize(mat, (len2+1));   // One column for the diagonal matrix
    resize(horizontal, (len2+1));   // One column for the horizontal matrix
    resize(trace, len1 * ((len2 >> 1) + (len2 & 1)), 0);
    TTraceValue tvMat= 0;

    // Record the max score
    TScoreValue score_max = 0;
    indexPair[0] = 0; indexPair[1] = 0;
    initialDir = Stop;

    // Classical DP
    TScoreValue max_val = 0;
    TScoreValue a = 0;
    TScoreValue b = 0;
    typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
    TTraceIter it = begin(trace, Standard() );
    TMatIter matIt = begin(mat, Standard() );
    TMatIter horiIt = begin(horizontal, Standard() );
    *matIt = 0;
    for(TSize row = 1; row <= len2; ++row) {
        *(++matIt) = 0;
        *(++horiIt) = scoreGapOpenHorizontal(sc, sequenceEntryForScore(sc, str1, 0),
                                             sequenceEntryForScore(sc, str2, row - 1)) -
                     scoreGapExtendHorizontal(sc, sequenceEntryForScore(sc, str1, 0),
                                              sequenceEntryForScore(sc, str2, row - 1));
    }
    for(TSize col = 1; col <= len1; ++col) {
        matIt = begin(mat, Standard() );
        horiIt = begin(horizontal, Standard() );
        TScoreValue diagValMat = *matIt;
        *matIt = 0;
        vert = scoreGapOpenVertical(sc, sequenceEntryForScore(sc, str1, col-1), sequenceEntryForScore(sc, str2, 0)) -
               scoreGapExtendVertical(sc, sequenceEntryForScore(sc, str1, col-1), sequenceEntryForScore(sc, str2, 0));
        TSize row = 1;
        while(row <= len2) {
            if (_isClumping(forbidden, row, col, len2)) {
                *it <<= 3;
                *it |= Stop;
                max_val = 0;
                vert = 0;
                *(++horiIt) = 0;
                ++matIt;
            } else {
                // Get the new maximum for vertical
                a = *matIt + scoreGapOpenVertical(sc, sequenceEntryForScore(sc, str1, col-1),
                                                  sequenceEntryForScore(sc, str2, row-1));
                b = vert + scoreGapExtendVertical(sc, sequenceEntryForScore(sc, str1, col-1),
                                                  sequenceEntryForScore(sc, str2, row-1));
                if (a > b) { vert = a; *it |= 1;}
                else vert = b;

                // Get the new maximum for horizontal
                *it <<= 1;
                a = *(++matIt) + scoreGapOpenHorizontal(sc, sequenceEntryForScore(sc, str1, col-1),
                                                        sequenceEntryForScore(sc, str2, row-1));
                b = *(++horiIt) + scoreGapExtendHorizontal(sc, sequenceEntryForScore(sc, str1, col-1),
                                                           sequenceEntryForScore(sc, str2, row-1));
                if (a > b) {*horiIt = a; *it |= 1; }
                else *horiIt =  b;

                // Get the new maximum for mat
                *it <<= 2;
                max_val = diagValMat + score(const_cast<TScore&>(sc), sequenceEntryForScore(sc, str1, col-1),
                                             sequenceEntryForScore(sc, str2, row-1));
                tvMat =  Diagonal;
                if (vert > max_val) {
                    max_val = vert;
                    tvMat =  Vertical;
                }
                if (*horiIt > max_val) {
                    max_val = *horiIt;
                    tvMat =  Horizontal;
                }
                if (0 >= max_val) {
                    max_val = 0;
                    tvMat =  Stop;
                }
                *it |= tvMat;
            }

            // Assign the new diagonal values
            diagValMat = *matIt;
            *matIt = max_val;

            // Record the new best score
            if (max_val > score_max) {
                indexPair[0] = row; indexPair[1] = col;
                score_max = max_val;
                initialDir = tvMat;
            }

            if (row & 1) *it <<= 1; else ++it;
            ++row;
        }
        if (!(row & 1)) {*it <<= 3; ++it; }
    }

    //// Debug code
    //for(TSize i= 0; i<len2;++i) {
    //    for(TSize j= 0; j<len1;++j) {
    //        std::cout << (TSize) getValue(trace, j*len2 + i) << ',';
    //    }
    //    std::cout << std::endl;
    //}
    //std::cout << "Max score: " << best_row << ',' << best_col << ':' << score_max << " (" << (TSize) initialDir << ")" << std::endl;

    return score_max;
}

template<typename TAlign, typename TStringSet, typename TForbidden, typename TScore>
inline typename Value<TScore>::Type
_localAlignment(TAlign& align,
                TStringSet& str,
                TForbidden& forbidden,
                TScore const& sc,
                SmithWatermanClump)
{
    SEQAN_CHECKPOINT
    typedef typename Value<TScore>::Type TScoreValue;
    typedef typename Size<TStringSet>::Type TSize;

    TScoreValue maxScore;
    TSize indexPair[2];

    // Trace
    String<unsigned char> trace;
    unsigned char initialDir;

    // Create the trace
    maxScore = _alignSmithWaterman(trace, str, sc, initialDir, indexPair, forbidden);

    //// Debug code
    //for(TSize i= 0; i<length(str[1]);++i) {
    //    for(TSize j= 0; j<length(str[0]);++j) {
    //        std::cout << (TSize) getValue(forbidden, j*length(str[1]) + i) << ',';
    //    }
    //    std::cout << std::endl;
    //}
    //std::cout << std::endl;

    // Follow the trace and create the alignment
    _alignSmithWatermanTrace(align, str, trace, initialDir, indexPair, forbidden);

    return maxScore;
}

template<typename TString, typename TMatches, typename TScores, typename TScore, typename TSize1>
inline void
_localAlignment(StringSet<TString, Dependent<> > const& str,
                TMatches& matches,
                TScores& scores,
                TScore const& sc,
                TSize1 numAlignments,
                SmithWatermanClump)
{
    SEQAN_CHECKPOINT
    typedef typename Value<TScore>::Type TScoreValue;
    typedef typename Size<TMatches>::Type TSize;

    // For clumpping remember the used positions
    TSize len0 = length(str[0]);
    TSize len1 = length(str[1]);
    String<bool> forbidden;
    resize(forbidden, len0 * len1, false);

    // Stop looking for local alignments, if there score is too low
    TScoreValue local_score = 0;
    TScoreValue last_score = 0;
    for(TSize count = 0; count < (TSize) numAlignments; ++count) {
        // Create the local alignment
        TSize from = length(matches);
        local_score = _localAlignment(matches, str, forbidden, sc, SmithWatermanClump());
        TSize to = length(matches);
        if (2 * local_score < last_score) {
            resize(matches, from, Generous());
            break;
        }
        last_score = local_score;
        resize(scores, to);
        for(TSize k = from; k<to; ++k) scores[k] = local_score;
    }
}

template<typename TString, typename TMatches, typename TScores, typename TScoreValue, typename TSpec2, typename TSize, typename TTag>
inline void
_multiLocalAlignment(StringSet<TString, Dependent<> > const& str,
                     TMatches& matches,
                     TScores& scores,
                     Score<TScoreValue, TSpec2> const& sc,
                     TSize numAlignments,
                     TTag)
{
    // Make a multiple local alignment and get all matches
    _localAlignment(str,matches,scores,sc,numAlignments,TTag());
}

// ---------------------------------------------------------------------------
// END OF TO BE REMOVED LEGACY CODE
// ---------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////////////
// Pair selection to calculate alignment
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

// Dummy function selecting all pairs
template<typename TString, typename TSpec, typename TSize2, typename TSpec2>
inline void
selectPairs(StringSet<TString, TSpec> const& str,
            String<TSize2, TSpec2>& pList)
{
    SEQAN_CHECKPOINT
    typedef StringSet<TString, TSpec> TStringSet;
    typedef typename Size<TStringSet>::Type TSize;
    typedef typename Iterator<String<TSize2, TSpec2>, Standard>::Type TPairIter;

    TSize nseq = length(str);
    resize(pList, nseq * (nseq - 1));
    TPairIter itPair = begin(pList, Standard());
    for(TSize i=0; i<nseq-1; ++i) {
        for(TSize j=i+1; j<nseq; ++j) {
            *itPair = i; ++itPair;
            *itPair = j; ++itPair;
        }
    }
}


//////////////////////////////////////////////////////////////////////////////
// Alignment statistics
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

// TODO(holtgrew): Hard-code TSize1 as __int64, size_t?

template<typename TFragment, typename TSpec1, typename TStringSet, typename TPos, typename TSize1>
inline void
getAlignmentStatistics(String<TFragment, TSpec1> const& matches,
                       TStringSet& str,
                       TPos const from,
                       TPos const to,
                       TSize1& matchLength,    // Number of identical characters
                       TSize1& overlapLength,    // Number of character in overlapping segments (with mismatches and gaps)
                       TSize1& alignLength)    // Length of the alignment
{
    typedef String<TFragment, TSpec1> const TFragmentMatches;
    typedef typename Size<TFragmentMatches>::Type TSize;
    typedef typename Id<TFragmentMatches>::Type TId;
    typedef typename Iterator<TFragmentMatches, Standard>::Type TFragIter;
    typedef typename Value<TStringSet>::Type TString;
    typedef typename Value<TString>::Type TAlphabet;
    matchLength = 0;
    TSize len1 = length(str[0]);
    TSize len2 = length(str[1]);


    TSize minId1 = len1 + len2;
    TSize minId2 = len1 + len2;
    TSize maxId1 = 0;
    TSize maxId2 = 0;
    TSize matchMismatch_length = 0;

    TFragIter itFrag = begin(matches, Standard());
    TFragIter itFragEnd = itFrag;
    itFrag += from;
    itFragEnd += to;
    TId id1 = sequenceId(*itFrag, 0);
    TId id2 = sequenceId(*itFrag, 1);
    TSize fragLen = 0;
    TSize beginI = 0;
    TSize beginJ = 0;
    for(;itFrag != itFragEnd; ++itFrag) {
        fragLen = fragmentLength(*itFrag, id1);
        beginI = fragmentBegin(*itFrag, id1);
        beginJ = fragmentBegin(*itFrag, id2);
        if (beginI < minId1) minId1 = beginI;
        if (beginJ < minId2) minId2 = beginJ;
        if (beginI + fragLen > maxId1) maxId1 = beginI + fragLen;
        if (beginJ + fragLen > maxId2) maxId2 = beginJ + fragLen;
        typedef typename Infix<TString>::Type TInfix;
        typedef typename Iterator<TInfix, Standard>::Type TInfixIter;
        TInfix inf1 = label(*itFrag, str, id1);
        TInfix inf2 = label(*itFrag, str, id2);
        TInfixIter sIt1 = begin(inf1, Standard());
        TInfixIter sIt2 = begin(inf2, Standard());
        TInfixIter sIt1End = end(inf1, Standard());
        matchMismatch_length += fragLen;
        for(;sIt1 != sIt1End; ++sIt1, ++sIt2)
            if ( (TAlphabet) *sIt1  == (TAlphabet) *sIt2) ++matchLength;
    }
    alignLength = static_cast<TSize1>(matchMismatch_length + (len1 - matchMismatch_length) + (len2 - matchMismatch_length));
    overlapLength = alignLength -  minId1 - minId2 - (len1 + len2 - maxId1 - maxId2);
}

//////////////////////////////////////////////////////////////////////////////
// Segment Match Generation
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TSize2, typename TSpec2, typename TSegmentMatches, typename TScores>
inline void
appendSegmentMatches(StringSet<TString, Dependent<TSpec> > const& str,
                     String<TSize2, TSpec2> const& pList,
                     TSegmentMatches& matches,
                     TScores& scores,
                     LcsLibrary)
{
    typedef StringSet<TString, Dependent<TSpec> > TStringSet;
    typedef String<TSize2, TSpec2> TPairList;
    typedef typename Size<TStringSet>::Type TSize;
    typedef typename Id<TStringSet>::Type TId;
    //typedef typename Value<TSegmentMatches>::Type TFragment;
    //typedef typename Value<TScores>::Type TScoreValue;
    typedef typename Iterator<TPairList const, Standard>::Type TPairIter;

    // Pairwise longest common subsequence
    TPairIter itPair = begin(pList, Standard());
    TPairIter itPairEnd = end(pList, Standard());
    for(;itPair != itPairEnd; ++itPair) {
        TStringSet pairSet;
        TId id1 = positionToId(str, *itPair); ++itPair;
        TId id2 = positionToId(str, *itPair);
        assignValueById(pairSet, const_cast<TStringSet&>(str), id1);
        assignValueById(pairSet, const_cast<TStringSet&>(str), id2);

        // Lcs between first and second string
        TSize from = length(matches);
        globalAlignment(matches, pairSet, Lcs());
        TSize to = length(matches);

        // Record the scores
        resize(scores, to);
        typedef typename Iterator<TSegmentMatches, Standard>::Type TMatchIter;
        typedef typename Iterator<TScores, Standard>::Type TScoreIter;
        TScoreIter itScore = begin(scores, Standard());
        TScoreIter itScoreEnd = end(scores, Standard());
        TMatchIter itMatch = begin(matches, Standard());
        itScore+=from;
        itMatch+=from;
        for(;itScore != itScoreEnd; ++itScore, ++itMatch) *itScore = (*itMatch).len;
    }
}


//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TSegmentMatches, typename TScores, typename TAlphabet, typename TSize>
inline void
appendSegmentMatches(StringSet<TString, Dependent<TSpec> > const& str,
                     TSegmentMatches& matches,
                     TScores& scores,
                     TSize ktup,
                     TAlphabet,
                     KmerLibrary)
{
    //typedef StringSet<TString, Dependent<TSpec> > TStringSet;
    typedef typename Value<TScores>::Type TScoreValue;
    typedef typename Value<TSegmentMatches>::Type TFragment;
    //typedef typename Id<TStringSet>::Type TId;
    typedef String<TSize> TTupelString;
    typedef String<TTupelString> TTupelStringSet;

    // Initialization
    TSize nseq = length(str);
    TSize alphabet_size = ValueSize<TAlphabet>::VALUE;

    // Transform the set of strings into a set of strings of k-tupels
    TTupelStringSet tupSet;
    resize(tupSet, nseq);
    for(TSize k=0;k<nseq;++k) {
        _getTupelString(str[k], tupSet[k], ktup, TAlphabet());
    }

    // Build one q-gram Index for all sequences
    typedef std::pair<TSize, TSize> TPosSeqPair;
    typedef std::set<TPosSeqPair> TQGramOcc;
    String<TQGramOcc> qIndex;
    TSize qIndexSize = 1;
    for(TSize i=0; i<(TSize) ktup;++i) qIndexSize *= alphabet_size;
    resize(qIndex, qIndexSize);
    for(TSize k=0;k<nseq;++k) {
        for(TSize i = 0;i < (TSize) length(tupSet[k]);++i) {
            qIndex[ tupSet[k][i] ].insert(std::make_pair(i, k));
        }
    }
    for(TSize q=0;q< (TSize) qIndexSize;++q) {
        typename TQGramOcc::const_iterator pos = qIndex[q].begin();
        typename TQGramOcc::const_iterator posEnd = qIndex[q].end();
        while (pos != posEnd) {
            typename TQGramOcc::const_iterator pos2 = pos;
            ++pos2;
            while (pos2 != posEnd) {
                if (pos->second != pos2->second) {
                    appendValue(matches, TFragment(positionToId(str, pos->second), pos->first, positionToId(str, pos2->second), pos2->first, ktup));
                    appendValue(scores, (TScoreValue) ktup);
                }
                ++pos2;
            }
            ++pos;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TSegmentMatches, typename TScores, typename TSize>
inline void
appendSegmentMatches(StringSet<TString, Dependent<TSpec> > const& str,
                     TSegmentMatches& matches,
                     TScores& scores,
                     TSize ktup,
                     KmerLibrary)
{
    SEQAN_CHECKPOINT
    appendSegmentMatches(str, matches, scores, ktup,  typename Value<TString>::Type(), KmerLibrary());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TSegmentMatches, typename TScores>
inline void
appendSegmentMatches(StringSet<TString, Dependent<TSpec> > const& str,
                     TSegmentMatches& matches,
                     TScores& scores,
                     KmerLibrary)
{
    SEQAN_CHECKPOINT
    appendSegmentMatches(str, matches, scores, 3, KmerLibrary());
}


//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TSize2, typename TSpec2, typename TScore, typename TSegmentMatches, typename TScores>
inline void
appendSegmentMatches(StringSet<TString, Dependent<TSpec> > const& str,
                     String<TSize2, TSpec2> const& pList,
                     TScore const& score_type,
                     TSegmentMatches& matches,
                     TScores& scores,
                     LocalPairwiseLibrary)
{
    typedef StringSet<TString, Dependent<TSpec> > TStringSet;
    typedef String<TSize2, TSpec2> TPairList;
    //typedef typename Size<TStringSet>::Type TSize;
    typedef typename Id<TStringSet>::Type TId;
    typedef typename Iterator<TPairList const, Standard>::Type TPairIter;

    // Pairwise alignments
    TPairIter itPair = begin(pList, Standard());
    TPairIter itPairEnd = end(pList, Standard());
    for(;itPair != itPairEnd; ++itPair) {
        // Make a pairwise string-set
        TStringSet pairSet;
        TId id1 = positionToId(str, *itPair); ++itPair;
        TId id2 = positionToId(str, *itPair);
        assignValueById(pairSet, const_cast<TStringSet&>(str), id1);
        assignValueById(pairSet, const_cast<TStringSet&>(str), id2);

        _multiLocalAlignment(pairSet, matches, scores, score_type, 4, SmithWatermanClump());
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSpec, typename TSize>
inline void
_resizeWithRespectToDistance(String<TValue, TSpec>& dist,
                              TSize nseq)
{
    SEQAN_CHECKPOINT;
    resize(dist, nseq * nseq, 0);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TSize>
inline void
_resizeWithRespectToDistance(Graph<Undirected<TCargo, TSpec> >& dist, TSize nseq)
{
    SEQAN_CHECKPOINT
    clear(dist);
    reserve(_getVertexString(dist), nseq);
    for(TSize i=0;i<nseq; ++i) addVertex(dist);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSize>
inline void
_resizeWithRespectToDistance(Nothing&, TSize)
{
    SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFragment, typename TSpec1, typename TString, typename TSpec2, typename TValue,  typename TSpec, typename TSize>
inline void
_setDistanceValue(String<TFragment, TSpec1>& matches,
                   StringSet<TString, TSpec2>& pairSet,
                   String<TValue, TSpec>& dist,
                   TSize i,
                   TSize j,
                   TSize nseq,
                   TSize from)
{
    SEQAN_CHECKPOINT
    typedef typename Position<String<TFragment, TSpec1> >::Type TPos;

    // Determine a sequence weight
    TValue matchLen = 0;
    TValue overlapLen = 0;
    TValue alignLen = 0;
    getAlignmentStatistics(matches, pairSet, (TPos) from, (TPos) length(matches),  matchLen, overlapLen, alignLen);

    // Calculate sequence similarity
    TValue x = SEQAN_DISTANCE_UNITY - static_cast<TValue>(static_cast<double>(matchLen) / static_cast<double>(alignLen) * static_cast<double>(SEQAN_DISTANCE_UNITY));
    if (i < j)
      dist[i * nseq + j] = x;
    else
      dist[j * nseq + i] = x;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFragment, typename TSpec1, typename TString, typename TSpec2, typename TCargo,  typename TSpec, typename TSize>
inline void
_setDistanceValue(String<TFragment, TSpec1>& matches,
                   StringSet<TString, TSpec2>& pairSet,
                   Graph<Undirected<TCargo, TSpec> >& dist,
                   TSize i,
                   TSize j,
                   TSize,
                   TSize from)
{
    SEQAN_CHECKPOINT

    // Determine a sequence weight
    TCargo matchLen = 0;
    TCargo overlapLen = 0;
    TCargo alignLen = 0;
    getAlignmentStatistics(matches, pairSet, (TSize) from, (TSize) length(matches),  matchLen, overlapLen, alignLen);

    // Calculate sequence similarity
    TCargo normalizedSimilarity = SEQAN_DISTANCE_UNITY - (TCargo) (((double) matchLen / (double) overlapLen) * ((double) overlapLen / (double) alignLen) * (double) SEQAN_DISTANCE_UNITY);

    addEdge(dist, i, j, normalizedSimilarity);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TFragment, typename TSpec, typename TString, typename TSpec2, typename TSize>
inline void
_setDistanceValue(String<TFragment, TSpec>&,
                   StringSet<TString, TSpec2>&,
                   Nothing&,
                   TSize,
                   TSize,
                   TSize,
                   TSize)
{
    SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TSize2, typename TSpec2, typename TScore, typename TSegmentMatches, typename TScoreValues, typename TDistance, typename TAlignConfig>
inline void
appendSegmentMatches(StringSet<TString, Dependent<TSpec> > const& str,
                     String<TSize2, TSpec2> const& pList,
                     TScore const& score_type,
                     TSegmentMatches& matches,
                     TScoreValues& scores,
                     TDistance& dist,
                     TAlignConfig const& ac,
                     GlobalPairwiseLibrary)
{
    SEQAN_CHECKPOINT
    typedef StringSet<TString, Dependent<TSpec> > TStringSet;
    typedef typename Id<TStringSet>::Type TId;
    typedef typename Size<TStringSet>::Type TSize;
    typedef typename Value<TScoreValues>::Type TScoreValue;
    typedef typename Iterator<String<TSize2, TSpec2> const, Standard>::Type TPairIter;

    // Initialization
    TSize nseq = length(str);
    _resizeWithRespectToDistance(dist, nseq);

    // Pairwise alignments
    TPairIter itPair = begin(pList, Standard());
    TPairIter itPairEnd = end(pList, Standard());
    for(;itPair != itPairEnd; ++itPair) {
        // Make a pairwise string-set
        TStringSet pairSet;
        TId id1 = positionToId(str, *itPair); ++itPair;
        TId id2 = positionToId(str, *itPair);
        assignValueById(pairSet, const_cast<TStringSet&>(str), id1);
        assignValueById(pairSet, const_cast<TStringSet&>(str), id2);

        // Alignment
        TSize from = length(matches);
        TScoreValue myScore = globalAlignment(matches, pairSet, score_type, ac, Gotoh() );
        TSize to = length(matches);

        // Record the scores
        resize(scores, to);
        typedef typename Iterator<TScoreValues, Standard>::Type TScoreIter;
        TScoreIter itScore = begin(scores, Standard());
        TScoreIter itScoreEnd = end(scores, Standard());
        itScore+=from;
        for(;itScore != itScoreEnd; ++itScore) *itScore = myScore;

        // Get the alignment statistics
        _setDistanceValue(matches, pairSet, dist, (TSize) *(itPair-1), (TSize) *itPair, (TSize) nseq, (TSize)from);
    }
}


//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TSize2, typename TSpec2, typename TScore, typename TSegmentMatches, typename TScoreValues, typename TDistance>
inline void
appendSegmentMatches(StringSet<TString, Dependent<TSpec> > const& str,
                     String<TSize2, TSpec2> const& pList,
                     TScore const& score_type,
                     TSegmentMatches& matches,
                     TScoreValues& scores,
                     TDistance& dist,
                     GlobalPairwiseLibrary)
{
    SEQAN_CHECKPOINT
    appendSegmentMatches(str, pList, score_type, matches, scores, dist, AlignConfig<>(), GlobalPairwiseLibrary() );
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
