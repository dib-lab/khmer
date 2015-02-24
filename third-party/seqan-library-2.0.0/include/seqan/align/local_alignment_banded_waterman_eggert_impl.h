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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// Author: Anne-Katrin Emde <anne-katrin.emde@fu-berlin.de>
// Author: Tobias Rausch <rausch@embl.de>
// Author: Birte Kehr <birte.kehr@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_BANDED_WATERMAN_EGGERT_IMPL_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_BANDED_WATERMAN_EGGERT_IMPL_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Legacy.
struct BandedWatermanEggert_;
typedef Tag<BandedWatermanEggert_> BandedWatermanEggert;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _alignBandedSmithWaterman()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceH, typename TSequenceV, typename TScore, typename TDiagonal>
inline TScoreValue
_alignBandedSmithWaterman(LocalAlignmentFinder<TScoreValue>& finder,
                          TSequenceH const & seqH,
                          TSequenceV const & seqV,
                          TScore const& sc,
                          TScoreValue const cutoff,
                          TDiagonal const diagL,
                          TDiagonal const diagU) {
    typedef typename Size<TSequenceH>::Type TSize;
    typedef LocalAlignmentFinder<TScoreValue> TFinder;

    // Initialization
    TSequenceH const& str1 = seqH;
    TSequenceV const& str2 = seqV;
    TSize len1 = length(str1);
    TSize len2 = length(str2);

    TSize lo_row = (diagU <= 0) ? static_cast<TSize>(-diagU) : 0;
    TSize hi_row = len2;// + 1;
    if (len1 - diagL < hi_row) hi_row = static_cast<TSize>(len1 - diagL);

    TSize height = hi_row - lo_row + 1;
    TSize diagonalWidth = (TSize) (diagU - diagL + 1);

    TSize actualCol, actualRow;
    TScoreValue verti_val, hori_val;

    // Initialize iterators
    typename TFinder::TMatrixIterator matIt = begin(finder.matrix);  // Iterator in current row
    goNext(matIt, 1);
    typename TFinder::TMatrixIterator matIt2 = begin(finder.matrix); // Iterator in previous row (for diagonal and vertical value)

    for (TSize row = 1; row < height; ++row) {
        actualRow = row + lo_row;
        hori_val = 0;

        for (TSize col = 0; col < diagonalWidth; ++col, ++matIt) {
            // handle begin and end triangle of band
            if ((int)col + diagL + (int)actualRow < 0) {++matIt2; continue;}
            actualCol = static_cast<TSize>(col + diagL + actualRow);
            if (actualCol > len1) {++matIt2; continue;}

            if (actualCol != 0) {
                // Get the new maximum for diagonal
                *matIt = *matIt2 +
                            score(const_cast<TScore&>(sc), sequenceEntryForScore(sc, str1, (int) actualCol - 1),
                                  sequenceEntryForScore(sc, str2, (int) actualRow - 1));

                ++matIt2;

                // Get the new maximum for vertical
                if (col < diagonalWidth - 1) {
                    verti_val = *matIt2 +
                                    scoreGapExtendVertical(sc, sequenceEntryForScore(sc, str1, (int) actualCol - 1),
                                                           sequenceEntryForScore(sc, str2, (int) actualRow - 1));
                    if (verti_val > *matIt) {
                        *matIt = verti_val;
                    }
                }

                // Get the new maximum for horizontal
                if (col > 0) {
                    hori_val = hori_val +
                                scoreGapExtendHorizontal(sc, sequenceEntryForScore(sc, str1, (int) actualCol - 1),
                                                         sequenceEntryForScore(sc, str2, (int) actualRow - 1));
                    if (hori_val > *matIt) {
                        *matIt = hori_val;
                    }
                }

                // Check if new maximum is greater than 0
                if (0 > *matIt) {
                    *matIt = 0;
                }

                // Record the new best score
                if (*matIt >= cutoff) {
                    push(finder.pQ, ScoreAndID<TScoreValue, typename TFinder::TMatrixPosition>(*matIt, position(matIt)));
                }
            } else {
                // First column (*matIt = 0)
                ++matIt2;
            }
            hori_val = *matIt;
        }
    }
 //   // Debug code
 //   std::cerr << std::endl;
    //for(TSize i= 0; i<height; ++i) {
    //    for(TSize j= 0; j<diagonalWidth; ++j) {
    //        std::cerr << value(finder.matrix, j, i) << ',';
    //    }
    //    if (i > 0) std::cerr << " " << str2[i-1] << std::endl;
    //    else std::cerr << std::endl;
    //}
 //   if(length(finder.pQ) > 0) std::cerr << "Max score: " << top(finder.pQ).value_ << std::endl;

    if(!empty(finder.pQ)) {
        finder.bestEndPos = top(finder.pQ).id_;
        return top(finder.pQ).value_;
    } else {
        return 0;
    }
}

// ----------------------------------------------------------------------------
// Helper Function _alignBandedSmithWatermanDeclump()
// ----------------------------------------------------------------------------

template<typename TSpec, typename TSize>
inline void
_setForbiddenCell(String<bool, TSpec>& forbidden,
                  TSize len1,
                  TSize len2,
                  TSize numRows)
{
    forbidden[(len1 - 1)*numRows + (len2 - 1)] = true;
}


template<typename TSize>
inline void
_setForbiddenCell(Nothing&,
                  TSize,
                  TSize,
                  TSize)
{
}

// ----------------------------------------------------------------------------
// Function _alignBandedSmithWatermanDeclump()
// ----------------------------------------------------------------------------

template<typename TScoreValue, typename TSequenceH, typename TSequenceV, typename TScore, typename TDiagonal>
inline void
_alignBandedSmithWatermanDeclump(LocalAlignmentFinder<TScoreValue>& finder,
                                 TSequenceH const & seqH,
                                 TSequenceV const & seqV,
                                 TScore const& sc,
                                 TScoreValue const cutoff,
                                 TDiagonal const diagL,
                                 TDiagonal const diagU)
{
    typedef typename Size<TSequenceH>::Type TSize;
    typedef LocalAlignmentFinder<TScoreValue> TFinder;
    typedef unsigned char TTraceValue;

    // Traceback values
    TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;

    TSequenceH const& str1 = seqH;
    TSequenceV const& str2 = seqV;
    TSize len1 = length(str1);

    TSize diagonalWidth = (TSize) (diagU - diagL + 1);
    TSize lo_row = (diagU <= 0) ? static_cast<TSize>(-diagU) : 0;
    TSize hi_row = length(seqV);// + 1;
    if (len1 - diagL < hi_row) hi_row = static_cast<TSize>(len1 - diagL);
    TSize height = hi_row - lo_row + 1;

    TSize actualRow, actualCol;
    //TSize row, col;

    // Initialize iterators
    typename TFinder::TMatrixIterator matIt = iter(finder.matrix, finder.bestBeginPos);
    typename TFinder::TMatrixIterator matIt2;

    // Initialize column boundaries
    TSize minCol = diagonalWidth, newMinCol = diagonalWidth;
    TSize maxCol = 0, newMaxCol = 0;
    TSize row = coordinate(matIt, 1);
    TSize col = coordinate(matIt, 0);
    TSize traceCol = col;

    matIt -= col;

    // Initialize position in trace
    TSize tracePos = length(finder.trace.sizes);
    while (tracePos > 0 && finder.trace.tvs[tracePos-1] != Diagonal) --tracePos;
    if (tracePos == 0) return;
    TSize traceSize = finder.trace.sizes[tracePos-1];
    TTraceValue traceValue = finder.trace.tvs[tracePos-1];

    // iterate over rows
    while((row < height) && ((maxCol > minCol) || (tracePos > 0))) {
        actualRow = row + lo_row;

        // make sure that all matrix entries of trace are re-calculated and set to forbidden
        while (traceSize == 0 && tracePos > 1) {
            // determine next trace direction
            --tracePos;
            traceValue = finder.trace.tvs[tracePos-1];
            if (traceValue == Horizontal) {
                traceCol += finder.trace.sizes[tracePos-1];
            } else {
                traceSize = finder.trace.sizes[tracePos-1];
            }
        }
        if (tracePos > 0) {
            // follow the trace in the current row
            if (traceValue == Diagonal) {
                --traceSize;
                _setForbiddenCell(finder.forbidden, row+1, traceCol+1, diagonalWidth);
                minCol = _min(minCol, traceCol);
                maxCol = _max(maxCol, traceCol+1);
            } else if (traceValue == Vertical) {
                if (traceCol > 0) --traceCol;
                --traceSize;
            }

            if (traceCol >= maxCol) {
                    maxCol = traceCol + 1;
            }
        }

        // iterate over columns that have to be re-calculated
        if (maxCol > minCol) {
            col = minCol;
            matIt += col;
            matIt2 = matIt - diagonalWidth;
            while (col < maxCol) {
                actualCol = static_cast<TSize>(col + diagL + actualRow);
                if (actualCol > len1) break;

                TScoreValue newVal = 0;

                // diagonal
                if (!value(finder.forbidden, position(matIt))) {
                    newVal = _max(newVal, *matIt2 + score(sc, sequenceEntryForScore(sc, str1, (int)actualCol-1),
                                                          sequenceEntryForScore(sc, str2, (int)actualRow-1)));
                }
                ++matIt2;

                // horizontal
                if (col > 0) {
                    newVal = _max(newVal, *(matIt-1) +
                                  scoreGapExtendHorizontal(sc, sequenceEntryForScore(sc, str1, (int) actualCol-1),
                                                           sequenceEntryForScore(sc, str2, (int)actualRow-1)));
                }

                // vertical
                if (col+1 < maxCol) {
                    newVal = _max(newVal, *matIt2 +
                                  scoreGapExtendVertical(sc, sequenceEntryForScore(sc, str1, (int)actualCol-1),
                                                         sequenceEntryForScore(sc, str2, (int)actualRow-1)));
                }

                if (newVal != *matIt) {
                    // matrix entry changed
                    *matIt = newVal;
                    maxCol = _min(_max(maxCol, col+2), diagonalWidth);
                    newMaxCol = _max(newMaxCol, col+1);
                    newMinCol = ((col != 0) ? _min(newMinCol, col-1) : 0);

                    // Record the new best score
                    if (newVal >= cutoff) {
                        push(finder.pQ, ScoreAndID<TScoreValue, typename TFinder::TMatrixPosition>(newVal, position(matIt)));
                    }
                } else {
                    // matrix entry did not change
                    if (col == minCol) {
                        ++minCol;
                    }
                }
                ++col;
                ++matIt;
            }
            matIt += diagonalWidth - col;
        } else {
            matIt += diagonalWidth;
        }

        minCol = _max((TSize)0, newMinCol);
        maxCol = _min(diagonalWidth, newMaxCol);
        newMinCol = diagonalWidth;
        newMaxCol = 0;

        ++row;
    }
  //  // Debug code
  //  std::cerr << std::endl;
  //  for (TSize i = 0; i < height; ++i) {
  //      for(TSize j = 0; j < diagonalWidth; ++j) {
  //          std::cerr << value(finder.matrix, j, i) << ',';
  //      }
        //if (i > 0) std::cerr << " " << str2[i-1] << "    ";
        //else std::cerr << "      ";
  //      for (TSize j= 0; j<diagonalWidth; ++j) {
  //          std::cerr << value(finder.forbidden, j+i*diagonalWidth) << ',';
  //      }
  //      if (i > 0) std::cerr << " " << str2[i-1] << std::endl;
        //else std::cerr << std::endl;
  //  }
}

// ----------------------------------------------------------------------------
// Function _alignBandedSmithWatermanTrace()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceH, typename TSequenceV, typename TId, typename TScore, typename TDiagonal>
inline Pair<Pair<TDiagonal> >
_alignBandedSmithWatermanTrace(LocalAlignmentFinder<TScoreValue> & finder,
                               TSequenceH const & seqH,
                               TSequenceV const & seqV,
                               TId id1,
                               TId id2,
                               TScore& sc,
                               TDiagonal diagL,
                               TDiagonal diagU)
{
    typedef LocalAlignmentFinder<TScoreValue> TFinder;
    typedef typename Size<TSequenceH>::Type TSize;
    typedef unsigned char TTraceValue;

    clear(finder.trace.sizes);
    clear(finder.trace.tvs);

    // Traceback values
    TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2; TTraceValue Stop = 3;

    // Initialization
    TSequenceH const& str1 = seqH;
    TSequenceV const& str2 = seqV;
    TSize len1 = length(str1);
    TSize len2 = length(str2);

    TSize lo_row = (diagU <= 0) ? static_cast<TSize>(-diagU) : 0;
    //TSize diagonalWidth = static_cast<TSize>(diagU - diagL + 1);

    // Start the trace from the cell with the max value
    typename TFinder::TMatrixIterator matIt = iter(finder.matrix, finder.bestEndPos);
    typename TFinder::TMatrixIterator matIt2 = iter(finder.matrix, finder.bestEndPos);
    goPrevious(matIt2, 1);

    TSize row = coordinate(matIt, 1);
    TSize col = coordinate(matIt, 0);
    TSize endRow = row + lo_row;
    TSize endCol = static_cast<TSize>(col + diagL + endRow);

    TSize actualRow = row + lo_row;
    TSize actualCol = static_cast<TSize>(col + diagL + actualRow);
    if ((actualCol == 0) || (actualRow == 0))
        return Pair<Pair<TDiagonal> >();

    if (actualCol < len1) _alignTracePrint(finder.trace, seqH, seqV, id1, actualCol, id2, actualRow, len1 - actualCol, Horizontal);
    if (actualRow < len2) _alignTracePrint(finder.trace, seqH, seqV, id1, actualCol, id2, actualRow, len2 - actualRow, Vertical);

    TTraceValue traceValue = Stop;
    TTraceValue nextTraceValue = Horizontal;
    TSize segLen = 0;

    while (nextTraceValue != Stop) {
        traceValue = nextTraceValue;
        if (*matIt == 0) {
            nextTraceValue = Stop;
        }
        else if (*matIt == *matIt2 +
            score(const_cast<TScore&>(sc), sequenceEntryForScore(sc, str1, (int)actualCol-1),
                  sequenceEntryForScore(sc, str2, (int)actualRow-1)))
        {
            nextTraceValue = Diagonal;
            --actualRow; --actualCol;
            --row;
            goPrevious(matIt, 1);
            goPrevious(matIt2, 1);
        } else if (*matIt == *(matIt2+1) +
            scoreGapExtendVertical(sc, sequenceEntryForScore(sc, str1, (int)actualCol-1),
                                   sequenceEntryForScore(sc, str2, (int)actualRow-1)))
        {
            nextTraceValue = Vertical;
            --actualRow;
            --row; ++col;
            goPrevious(matIt, 1); goNext(matIt, 0);
            goPrevious(matIt2, 1); goNext(matIt2, 0);
        } else {
            SEQAN_ASSERT_EQ(*matIt, *(matIt-1) +
                            scoreGapExtendHorizontal(sc, sequenceEntryForScore(sc, str1, (int) actualCol-1),
                                                     sequenceEntryForScore(sc, str2, (int) actualRow-1)));
            nextTraceValue = Horizontal;
            --actualCol;
            --col;
            goPrevious(matIt, 0);
            goPrevious(matIt2, 0);
        }
        if (traceValue == nextTraceValue) {
            ++segLen;
        } else {
            _alignTracePrint(finder.trace, seqH, seqV, id1, actualCol, id2, actualRow, segLen, traceValue);
            segLen = 1;
        }
    }

    // Handle the remaining sequence
    if (actualCol != 0) _alignTracePrint(finder.trace, seqH, seqV, (TId) id1, (TSize) 0, (TId) 0, (TSize) 0, (TSize) actualCol, Horizontal);
    if (actualRow != 0) _alignTracePrint(finder.trace, seqH, seqV, (TId) 0, (TSize) 0, (TId) id2, (TSize) 0, (TSize) actualRow, Vertical);

    goNext(matIt, 1); // assumes that each trace ends with a diagonal
    finder.bestBeginPos = position(matIt);

    return Pair<Pair<TDiagonal> >(Pair<TDiagonal>(actualCol, endCol), Pair<TDiagonal>(actualRow, endRow));
}

// ----------------------------------------------------------------------------
// Function _initLocalAlignmentFinder()
// ----------------------------------------------------------------------------

template<typename TSequenceH, typename TSequenceV, typename TScoreValue, typename TDiagonal>
void
_initLocalAlignmentFinder(TSequenceH const & seqH,
                          TSequenceV const & seqV,
                          LocalAlignmentFinder<TScoreValue> & finder,
                          BandedWatermanEggert const &,
                          TDiagonal const lowerDiag,
                          TDiagonal const upperDiag) {
    typedef LocalAlignmentFinder<TScoreValue> TFinder;
    typedef typename TFinder::TMatrix TMatrix;
    typedef typename Size<TMatrix>::Type TSize;

    TSize lo_row = (upperDiag <= 0) ? static_cast<TSize>(-upperDiag) : 0;
    TSize hi_row = length(seqV) + 1;
    TSize len0 = length(seqH);
    if (len0 - lowerDiag < hi_row) hi_row = static_cast<TSize>(len0 - lowerDiag);
    TSize height = hi_row - lo_row + 1;
    SEQAN_ASSERT_GEQ(upperDiag, lowerDiag);
    TSize diagonalWidth = (TSize) (upperDiag - lowerDiag + 1);

    setDimension(finder.matrix, 2);
    setLength(finder.matrix, 0, diagonalWidth);
    setLength(finder.matrix, 1, height);
    resize(finder.matrix, 0);

    resize(finder.forbidden, height * diagonalWidth, false);

    finder.bestEndPos = minValue<typename TFinder::TMatrixPosition>();
    finder.bestBeginPos = minValue<typename TFinder::TMatrixPosition>();
}

// ----------------------------------------------------------------------------
// Function _localAlignment()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceH, typename TGapsSpecH, typename TSequenceV, typename TGapsSpecV, typename TScoreValue2, typename TScoreSpec, typename TDiagonal>
inline TScoreValue
_localAlignment(LocalAlignmentFinder<TScoreValue> & finder,
                Gaps<TSequenceH, TGapsSpecH> & gapsH,
                Gaps<TSequenceV, TGapsSpecV> & gapsV,
                Score<TScoreValue2, TScoreSpec> const& sc,
                TScoreValue cutoff,
                TDiagonal diag1,
                TDiagonal diag2,
                BandedWatermanEggert const /*algorithmTag*/)
{
    clearGaps(gapsH);
    clearClipping(gapsH);
    clearGaps(gapsV);
    clearClipping(gapsV);

    _initLocalAlignmentFinder(source(gapsH), source(gapsV), finder, BandedWatermanEggert(), diag1, diag2);
    finder.needReinit = false;

    // Fill the matrix
    TScoreValue maxScore = _alignBandedSmithWaterman(finder, source(gapsH), source(gapsV), sc, cutoff, diag1, diag2);
    if (maxScore < cutoff) return 0;

    // Follow the matrix back from max entry and create a trace path
    Pair<Pair<TDiagonal> > alignmentPositions = _alignBandedSmithWatermanTrace(finder, source(gapsH), source(gapsV), 0u, 0u, sc, diag1, diag2);

    // Fill the gaps following the trace path.
    _pumpTraceToGaps(gapsH, gapsV, finder.trace);

    // The following (using rightOfGaps = false) only works if the gap open cost is negative.  In this case, we can
    // assume that there are no leading or trailing gaps in either sequence.
    //
    // We need this assumption because the we get source positions through alignmentPositions and conversion from
    // souce to view coordinates gives us the end of a gap if we are at the beginning.
    //
    // TODO(holtgrew): This should be done in a more robust way.  It would probably be better if alignmentPositions gave us the lengths of the local alignments!
    setClippedEndPosition(gapsH, toViewPosition(gapsH, alignmentPositions.i1.i2, false));
    setClippedEndPosition(gapsV, toViewPosition(gapsV, alignmentPositions.i2.i2, false));
    setClippedBeginPosition(gapsH, toViewPosition(gapsH, alignmentPositions.i1.i1));
    setClippedBeginPosition(gapsV, toViewPosition(gapsV, alignmentPositions.i2.i1));

    pop(finder.pQ);

    return maxScore;
}

// ----------------------------------------------------------------------------
// Function _localAlignmentNext()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceH, typename TGapsSpecH, typename TSequenceV, typename TGapsSpecV,typename TScore, typename TDiagonal>
inline TScoreValue
_localAlignmentNext(LocalAlignmentFinder<TScoreValue> & finder,
                    Gaps<TSequenceH, TGapsSpecH> & gapsH,
                    Gaps<TSequenceV, TGapsSpecV> & gapsV,
                    TScore const& sc,
                    TScoreValue cutoff,
                    TDiagonal diag1,
                    TDiagonal diag2,
                    BandedWatermanEggert const & /*algorithmTag*/) {
    // Declump the matrix and find new maximum score
    _alignBandedSmithWatermanDeclump(finder, source(gapsH), source(gapsV), sc, cutoff, diag1, diag2);
    typename LocalAlignmentFinder<TScoreValue>::TMatrixPosition nextBestEnd;
    nextBestEnd = _getNextBestEndPosition(finder, cutoff);
    if(nextBestEnd==0)
        return 0;
    TScoreValue maxScore = getValue(finder.matrix, nextBestEnd);
    if(maxScore == 0) return 0;

    // Follow the trace matrix and create a trace path
    Pair<Pair<TDiagonal> > alignmentPositions = _alignBandedSmithWatermanTrace(finder, source(gapsH), source(gapsV), 0u, 0u, sc, diag1, diag2);

    // Fill the gaps following the trace path.
    _pumpTraceToGaps(gapsH, gapsV, finder.trace);

    // The following (using rightOfGaps = false) only works if the gap open cost is negative.  In this case, we can
    // assume that there are no leading or trailing gaps in either sequence.
    //
    // We need this assumption because the we get source positions through alignmentPositions and conversion from
    // souce to view coordinates gives us the end of a gap if we are at the beginning.
    //
    // TODO(holtgrew): This should be done in a more robust way.  It would probably be better if alignmentPositions gave us the lengths of the local alignments!
    setClippedEndPosition(gapsH, toViewPosition(gapsH, alignmentPositions.i1.i2, false));
    setClippedEndPosition(gapsV, toViewPosition(gapsV, alignmentPositions.i2.i2, false));
    setClippedBeginPosition(gapsH, toViewPosition(gapsH, alignmentPositions.i1.i1));
    setClippedBeginPosition(gapsV, toViewPosition(gapsV, alignmentPositions.i2.i1));

    return maxScore;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_BANDED_WATERMAN_EGGERT_IMPL_H_
