// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Implements the affine gap cost functions.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_FORMULA_AFFINE_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_FORMULA_AFFINE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _retrieveTraceAffine                        [RecursionDirectionAll]
// ----------------------------------------------------------------------------

template <typename TScoreValue>
inline TraceBitMap_::TTraceValue
_retrieveTraceAffine(TScoreValue const & globalMax,
                     TScoreValue const & diagScore,
                     TScoreValue const & horiScore,
                     TScoreValue const & horiOpenScore,
                     TScoreValue const & vertiScore,
                     TScoreValue const & vertiOpenScore,
                     RecursionDirectionAll const &)
{
    typename TraceBitMap_::TTraceValue traceValue(TraceBitMap_::NONE);
    _conditionalOrOnEquality(traceValue, globalMax, diagScore, TraceBitMap_::DIAGONAL);
    _conditionalOrOnEquality(traceValue, globalMax, horiScore, TraceBitMap_::HORIZONTAL);
    _conditionalOrOnEquality(traceValue, globalMax, horiOpenScore, TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX);
    _conditionalOrOnEquality(traceValue, globalMax, vertiScore, TraceBitMap_::VERTICAL);
    _conditionalOrOnEquality(traceValue, globalMax, vertiOpenScore, TraceBitMap_::MAX_FROM_VERTICAL_MATRIX);
    return traceValue;
}

// ----------------------------------------------------------------------------
// Function _retrieveTraceAffine              [RecursionDirectionUpperDiagonal]
// ----------------------------------------------------------------------------

template <typename TScoreValue>
inline TraceBitMap_::TTraceValue
_retrieveTraceAffine(TScoreValue const & globalMax,
                     TScoreValue const & diagScore,
                     TScoreValue const & horiScore,
                     TScoreValue const & horiOpenScore,
                     TScoreValue const &,
                     TScoreValue const &,
                     RecursionDirectionUpperDiagonal const &)
{
    typename TraceBitMap_::TTraceValue traceValue(TraceBitMap_::NONE);
    _conditionalOrOnEquality(traceValue, globalMax, diagScore, TraceBitMap_::DIAGONAL);
    _conditionalOrOnEquality(traceValue, globalMax, horiScore, TraceBitMap_::HORIZONTAL);
    _conditionalOrOnEquality(traceValue, globalMax, horiOpenScore, TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX);
    return traceValue;
}

// ----------------------------------------------------------------------------
// Function _retrieveTraceAffine              [RecursionDirectionLowerDiagonal]
// ----------------------------------------------------------------------------

template <typename TScoreValue>
inline TraceBitMap_::TTraceValue
_retrieveTraceAffine(TScoreValue const & globalMax,
                     TScoreValue const & diagScore,
                     TScoreValue const &,
                     TScoreValue const &,
                     TScoreValue const & vertiScore,
                     TScoreValue const & vertiOpenScore,
                     RecursionDirectionLowerDiagonal const &)
{
    typename TraceBitMap_::TTraceValue traceValue(TraceBitMap_::NONE);
    _conditionalOrOnEquality(traceValue, globalMax, diagScore, TraceBitMap_::DIAGONAL);
    _conditionalOrOnEquality(traceValue, globalMax, vertiScore, TraceBitMap_::VERTICAL);
    _conditionalOrOnEquality(traceValue, globalMax, vertiOpenScore, TraceBitMap_::MAX_FROM_VERTICAL_MATRIX);
    return traceValue;
}

// ----------------------------------------------------------------------------
// Function _retrieveTraceAffine                 [RecursionDirectionHorizontal]
// ----------------------------------------------------------------------------

template <typename TScoreValue>
inline TraceBitMap_::TTraceValue
_retrieveTraceAffine(TScoreValue const & globalMax,
                     TScoreValue const &,
                     TScoreValue const & horiScore,
                     TScoreValue const & horiOpenScore,
                     TScoreValue const &,
                     TScoreValue const &,
                     RecursionDirectionHorizontal const &)
{
    typename TraceBitMap_::TTraceValue traceValue(TraceBitMap_::NONE);
    _conditionalOrOnEquality(traceValue, globalMax, horiScore, TraceBitMap_::HORIZONTAL);
    _conditionalOrOnEquality(traceValue, globalMax, horiOpenScore, TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX);
    return traceValue;
}

// ----------------------------------------------------------------------------
// Function _retrieveTraceAffine                   [RecursionDirectionVertical]
// ----------------------------------------------------------------------------

template <typename TScoreValue>
inline TraceBitMap_::TTraceValue
_retrieveTraceAffine(TScoreValue const & globalMax,
                     TScoreValue const &,
                     TScoreValue const &,
                     TScoreValue const &,
                     TScoreValue const & vertiScore,
                     TScoreValue const & vertiOpenScore,
                     RecursionDirectionVertical const &)
{
    typename TraceBitMap_::TTraceValue traceValue(TraceBitMap_::NONE);
    _conditionalOrOnEquality(traceValue, globalMax, vertiScore, TraceBitMap_::VERTICAL);
    _conditionalOrOnEquality(traceValue, globalMax, vertiOpenScore, TraceBitMap_::MAX_FROM_VERTICAL_MATRIX);
    return traceValue;
}

// ----------------------------------------------------------------------------
// Function _internalComputeScore      [RecursionDirectionDiagonal, AffineGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueGap>
inline typename TraceBitMap_::TTraceValue
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL,
                      TTraceValueGap,
                      TracebackOff const &,
                      RecursionDirectionDiagonal const &)
{
    if(activeCell._score < rightCompare)
        activeCell._score = rightCompare;
    return TraceBitMap_::NONE;
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueGap, typename TGapsPlacement>
inline typename TraceBitMap_::TTraceValue
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL leftTrace,
                      TTraceValueGap gapTrace,
                      TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> > const &,
                      RecursionDirectionDiagonal const &)
{
    if(activeCell._score <= rightCompare)
    {
        activeCell._score = rightCompare;
        return TraceBitMap_::DIAGONAL | leftTrace;
    }
    return leftTrace | gapTrace;
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueGap, typename TGapsPlacement>
inline typename TraceBitMap_::TTraceValue
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL leftTrace,
                      TTraceValueGap gapTrace,
                      TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> >  const &,
                      RecursionDirectionDiagonal const &)
{
    if(activeCell._score < rightCompare)
    {
        activeCell._score = rightCompare;  // Maximal score comes from diagonal.
        return TraceBitMap_::DIAGONAL | leftTrace;  // Return trace for Diagonal.
    }
    if (activeCell._score == rightCompare)      // Maximal score comes from previous computed directions and diagonal.
        return leftTrace | TraceBitMap_::DIAGONAL | gapTrace;   // Return all directions inclusively the flag indicating max from gap.
    return leftTrace | gapTrace; // Maximum comes from gap. Return gap value inclusively the flag indicating max from gap.
}

// ----------------------------------------------------------------------------
// Function _internalComputeScore    [RecursionDirectionHorizontal, AffineGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR>
inline typename TraceBitMap_::TTraceValue
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL,
                      TTraceValueR,
                      TracebackOff const &,
                      RecursionDirectionHorizontal const &)
{
    if(activeCell._horizontalScore < rightCompare)
        activeCell._score = activeCell._horizontalScore = rightCompare;
    else
        activeCell._score = activeCell._horizontalScore;
    return TraceBitMap_::NONE;
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR, typename TGapsPlacement>
inline typename TraceBitMap_::TTraceValue
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL leftTrace,
                      TTraceValueR  rightTrace,
                      TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> >  const &,
                      RecursionDirectionHorizontal const &)
{
    if(activeCell._horizontalScore < rightCompare)
    {
        activeCell._score = activeCell._horizontalScore = rightCompare;
        return rightTrace;
    }
    activeCell._score = activeCell._horizontalScore;
    return leftTrace;
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR, typename TGapsPlacement>
inline typename TraceBitMap_::TTraceValue
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL leftTrace,
                      TTraceValueR rightTrace,
                      TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> >  const &,
                      RecursionDirectionHorizontal const &)
{
    if(activeCell._horizontalScore < rightCompare)
    {
        activeCell._score = activeCell._horizontalScore = rightCompare;
        return rightTrace;
    }
    activeCell._score = activeCell._horizontalScore;
    if (activeCell._horizontalScore == rightCompare)
        return leftTrace | rightTrace;
    return leftTrace;
}

// ----------------------------------------------------------------------------
// Function _internalComputeScore      [RecursionDirectionVertical, AffineGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR>
inline typename TraceBitMap_::TTraceValue
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL,
                      TTraceValueR,
                      TracebackOff const &,
                      RecursionDirectionVertical const &)
{
    if(activeCell._verticalScore < rightCompare)
        activeCell._score = activeCell._verticalScore = rightCompare;
    else
        activeCell._score = activeCell._verticalScore;
    return TraceBitMap_::NONE;
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR, typename TGapsPlacement>
inline typename TraceBitMap_::TTraceValue
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL leftTrace,
                      TTraceValueR rightTrace,
                      TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> >  const &,
                      RecursionDirectionVertical const &)
{
    if(activeCell._verticalScore < rightCompare)
    {
        activeCell._score = activeCell._verticalScore = rightCompare;
        return rightTrace;
    }
    activeCell._score = activeCell._verticalScore;
    return leftTrace;
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR, typename TGapsPlacement>
inline typename TraceBitMap_::TTraceValue
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL leftTrace,
                      TTraceValueR rightTrace,
                      TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> >  const &,
                      RecursionDirectionVertical const &)
{
    if(activeCell._verticalScore < rightCompare)
    {
        activeCell._score = activeCell._verticalScore = rightCompare;
        return rightTrace;
    }
    activeCell._score = activeCell._verticalScore;
    if (activeCell._verticalScore == rightCompare)
        return leftTrace | rightTrace;
    return leftTrace;
}

// ----------------------------------------------------------------------------
// Function _internalComputeScore          [Vertical vs Horizontal, AffineGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue>
inline typename TraceBitMap_::TTraceValue
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TracebackOff const &)
{
    if(activeCell._score < activeCell._horizontalScore)
        activeCell._score = activeCell._horizontalScore;
    return TraceBitMap_::NONE;
}

template <typename TScoreValue, typename TGapsPlacement>
inline typename TraceBitMap_::TTraceValue
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> >  const &)
{
    if(activeCell._score < activeCell._horizontalScore)
    {
        activeCell._score = activeCell._horizontalScore;
        return TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX;
//        return traceRight;
    }
    return TraceBitMap_::MAX_FROM_VERTICAL_MATRIX;
//    returntraceRight;
}

template <typename TScoreValue, typename TGapsPlacement>
inline typename TraceBitMap_::TTraceValue
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> >  const &)
{
    if(activeCell._score < activeCell._horizontalScore)
    {
        activeCell._score = activeCell._horizontalScore;
        return TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX;
//        return traceRight | (traceLeft & TraceBitMap_::VERTICAL_OPEN);
    }
//    traceGap = TraceBitMap_::MAX_FROM_VERTICAL_MATRIX;
    if (activeCell._score == activeCell._horizontalScore)
    {
        return TraceBitMap_::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX;
//        return traceRight | traceLeft;
    }
    return TraceBitMap_::MAX_FROM_VERTICAL_MATRIX;
//    return traceLeft | (traceRight & TraceBitMap_::HORIZONTAL_OPEN);
}

// ----------------------------------------------------------------------------
// Function _doComputeScore                 [RecursionAllDirection, AffineGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig>
inline typename TraceBitMap_::TTraceValue
_doComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                DPCell_<TScoreValue, AffineGaps> const & previousDiagonal,
                DPCell_<TScoreValue, AffineGaps> const & previousHorizontal,
                DPCell_<TScoreValue, AffineGaps> const & previousVertical,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionAll const &,
                DPProfile_<TAlgorithm, AffineGaps, TTracebackConfig> const &)
{
    typedef typename TraceBitMap_::TTraceValue TTraceValue;

    // Now we have to find a smart version to solve this problem. Which is not as easy I would think.

    activeCell._horizontalScore = _horizontalScoreOfCell(previousHorizontal) +
                                        scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal);
    TScoreValue tmpScore = _scoreOfCell(previousHorizontal) + scoreGapOpenHorizontal(scoringScheme, seqHVal, seqVVal);

    TTraceValue tvGap = _internalComputeScore(activeCell, tmpScore, TraceBitMap_::HORIZONTAL, TraceBitMap_::HORIZONTAL_OPEN, TTracebackConfig(), RecursionDirectionHorizontal());
//    activeCell._score = activeCell._horizontalScore; _max(activeCell._horizontal, tmpScore)

    // Now we can decide for the optimal score in horizontal score or not?
    activeCell._verticalScore = _verticalScoreOfCell(previousVertical) + scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal);
    tmpScore = _scoreOfCell(previousVertical) + scoreGapOpenVertical(scoringScheme, seqHVal, seqVVal);
    tvGap |= _internalComputeScore(activeCell, tmpScore, TraceBitMap_::VERTICAL, TraceBitMap_::VERTICAL_OPEN, TTracebackConfig(), RecursionDirectionVertical());

    // Finds the maximum between the vertical and the horizontal matrix. Stores the flag for coming from a potential direction.
    TTraceValue tvMax = _internalComputeScore(activeCell, TTracebackConfig());  // Stores from where the maximal score comes.
    tmpScore = _scoreOfCell(previousDiagonal) + score(scoringScheme, seqHVal, seqVVal);
    return _internalComputeScore(activeCell, tmpScore, tvGap, tvMax, TTracebackConfig(), RecursionDirectionDiagonal());

    // This should be the fastest version.

//    TScoreValue tmpScoreDiagonal = _scoreOfCell(previousDiagonal) + score(scoringScheme, seqHVal, seqVVal);
//    TScoreValue tmpGapExtendHorizontal = _horizontalScoreOfCell(previousHorizontal) +
//                                         scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal);
//    TScoreValue tmpGapOpenHorizontal = _scoreOfCell(previousHorizontal)
//                                       + scoreGapOpenHorizontal(scoringScheme, seqHVal, seqVVal);
//    TScoreValue tmpGapExtendVertical = _verticalScoreOfCell(previousVertical)
//                                       + scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal);
//    TScoreValue tmpGapOpenVertical = _scoreOfCell(previousVertical)
//                                     + scoreGapOpenVertical(scoringScheme, seqHVal, seqVVal);
//
//    activeCell._horizontalScore = _max(tmpGapExtendHorizontal, tmpGapOpenHorizontal);
//    activeCell._verticalScore = _max(tmpGapExtendVertical, tmpGapOpenVertical);
//    activeCell._score = _max(tmpScoreDiagonal, _max(activeCell._horizontalScore, activeCell._verticalScore));
//
//    if (!IsTracebackEnabled_<TTracebackConfig>::VALUE)
//        return TraceBitMap_::NONE;
//
//    TTraceValue traceGapOpen = TraceBitMap_::NONE;
//    _conditionalOrOnInequality(traceGapOpen, _horizontalScoreOfCell(activeCell), tmpGapExtendHorizontal, TraceBitMap_::HORIZONTAL_OPEN);
//    _conditionalOrOnInequality(traceGapOpen, _verticalScoreOfCell(activeCell), tmpGapExtendVertical, TraceBitMap_::VERTICAL_OPEN);
//    return traceGapOpen |
//           _retrieveTraceAffine(activeCell._score, tmpScoreDiagonal, tmpGapExtendHorizontal, tmpGapOpenHorizontal,
//                                tmpGapExtendVertical, tmpGapOpenVertical, RecursionDirectionAll());
}

// ----------------------------------------------------------------------------
// Function _doComputeScore       [RecursionUpperDiagonalDirection, AffineGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig>
inline typename TraceBitMap_::TTraceValue
_doComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                DPCell_<TScoreValue, AffineGaps> const & previousDiagonal,
                DPCell_<TScoreValue, AffineGaps> const & previousHorizontal,
                DPCell_<TScoreValue, AffineGaps> const & /*previousVertical*/,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionUpperDiagonal const &,
                DPProfile_<TAlgorithm, AffineGaps, TTracebackConfig> const &)
{
    typedef typename TraceBitMap_::TTraceValue TTraceValue;
    activeCell._horizontalScore = _horizontalScoreOfCell(previousHorizontal)
                                         + scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal);

    activeCell._verticalScore = DPCellDefaultInfinity<DPCell_<TScoreValue, AffineGaps> >::VALUE;
    TScoreValue tmpScore = _scoreOfCell(previousHorizontal) + scoreGapOpenHorizontal(scoringScheme, seqHVal, seqVVal);
    TTraceValue tv = _internalComputeScore(activeCell, tmpScore, TraceBitMap_::HORIZONTAL, TraceBitMap_::HORIZONTAL_OPEN, TTracebackConfig(), RecursionDirectionHorizontal());
    tmpScore = _scoreOfCell(previousDiagonal) + score(scoringScheme, seqHVal, seqVVal);
    return _internalComputeScore(activeCell, tmpScore, tv, TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX, TTracebackConfig(), RecursionDirectionDiagonal());

//    TScoreValue tmpScoreDiagonal = _scoreOfCell(previousDiagonal) + score(scoringScheme, seqHVal, seqVVal);
//    TScoreValue tmpGapExtendHorizontal = _horizontalScoreOfCell(previousHorizontal)
//                                         + scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal);
//    TScoreValue tmpGapOpenHorizontal = _scoreOfCell(previousHorizontal)
//                                       + scoreGapOpenHorizontal(scoringScheme, seqHVal, seqVVal);
//
//    activeCell._horizontalScore = _max(tmpGapExtendHorizontal, tmpGapOpenHorizontal);
//    activeCell._score = _max(tmpScoreDiagonal, activeCell._horizontalScore);
//
//    if (!IsTracebackEnabled_<TTracebackConfig>::VALUE)
//        return TraceBitMap_::NONE;
//
//    TTraceValue traceGapOpen = TraceBitMap_::NONE;
//    _conditionalOrOnInequality(traceGapOpen, _horizontalScoreOfCell(activeCell), tmpGapExtendHorizontal, TraceBitMap_::HORIZONTAL_OPEN);
//    return traceGapOpen |
//           _retrieveTraceAffine(activeCell._score, tmpScoreDiagonal, tmpGapExtendHorizontal, tmpGapOpenHorizontal,
//                                TScoreValue(), TScoreValue(), RecursionDirectionUpperDiagonal());
}

// ----------------------------------------------------------------------------
// Function _doComputeScore       [RecursionDirectionLowerDiagonal, AffineGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig>
inline typename TraceBitMap_::TTraceValue
_doComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                DPCell_<TScoreValue, AffineGaps> const & previousDiagonal,
                DPCell_<TScoreValue, AffineGaps> const & /*previousHorizontal*/,
                DPCell_<TScoreValue, AffineGaps> const & previousVertical,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionLowerDiagonal const &,
                DPProfile_<TAlgorithm, AffineGaps, TTracebackConfig> const &)
{
    typedef typename TraceBitMap_::TTraceValue TTraceValue;

    activeCell._verticalScore = _verticalScoreOfCell(previousVertical) +
                                        scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal);
    TScoreValue tmpScore = _scoreOfCell(previousVertical) + scoreGapOpenVertical(scoringScheme, seqHVal, seqVVal);

    activeCell._horizontalScore = DPCellDefaultInfinity<DPCell_<TScoreValue, AffineGaps> >::VALUE;
    // This computes the difference between the vertical extend and vertical open.
    TTraceValue tv = _internalComputeScore(activeCell, tmpScore, TraceBitMap_::VERTICAL, TraceBitMap_::VERTICAL_OPEN, TTracebackConfig(), RecursionDirectionVertical());

    // Up to here, activeCell stores the highest value of vertical or vertical open.
    tmpScore = _scoreOfCell(previousDiagonal) + score(scoringScheme, seqHVal, seqVVal);
    return _internalComputeScore(activeCell, tmpScore, tv, TraceBitMap_::MAX_FROM_VERTICAL_MATRIX, TTracebackConfig(), RecursionDirectionDiagonal());  // Now we have this problem. How do we determine if the max comes from the vertical distance.

//    activeCell._verticalScore = _max(tmpGapExtendVertical, tmpGapOpenVertical);
//    activeCell._score = _max(tmpScoreDiagonal, activeCell._verticalScore);
//
//    if (!IsTracebackEnabled_<TTracebackConfig>::VALUE)
//        return TraceBitMap_::NONE;
//
//    TTraceValue traceGapOpen = TraceBitMap_::NONE;
//    _conditionalOrOnInequality(traceGapOpen, _verticalScoreOfCell(activeCell), tmpGapExtendVertical, TraceBitMap_::VERTICAL_OPEN);
//    return traceGapOpen |
//           _retrieveTraceAffine(activeCell._score, tmpScoreDiagonal, TScoreValue(), TScoreValue(), tmpGapExtendVertical,
//                                tmpGapOpenVertical, RecursionDirectionLowerDiagonal());
}

// ----------------------------------------------------------------------------
// Function _doComputeScore                      [RecursionHorizontalDirection]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig>
inline typename TraceBitMap_::TTraceValue
_doComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                DPCell_<TScoreValue, AffineGaps> const & /*previousDiagonal*/,
                DPCell_<TScoreValue, AffineGaps> const & previousHorizontal,
                DPCell_<TScoreValue, AffineGaps> const & /*previousVertical*/,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionHorizontal const &,
                DPProfile_<TAlgorithm, AffineGaps, TTracebackConfig> const &)
{
    //typedef typename TraceBitMap_::TTraceValue TTraceValue;
    TScoreValue tmpGapOpenHorizontal = _scoreOfCell(previousHorizontal) + scoreGapOpenHorizontal(scoringScheme, seqHVal, seqVVal);
    activeCell._horizontalScore = _horizontalScoreOfCell(previousHorizontal) + scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal);

    activeCell._verticalScore = DPCellDefaultInfinity<DPCell_<TScoreValue, AffineGaps> >::VALUE;
    return _internalComputeScore(activeCell, tmpGapOpenHorizontal, TraceBitMap_::HORIZONTAL, TraceBitMap_::HORIZONTAL_OPEN, TTracebackConfig(), RecursionDirectionHorizontal()) | TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX;

//    typedef typename TraceBitMap_::TTraceValue TTraceValue;
//    TScoreValue tmpGapOpenHorizontal = _scoreOfCell(previousHorizontal) +
//                                       scoreGapOpenHorizontal(scoringScheme, seqHVal, seqVVal);
//    TScoreValue tmpGapExtendHorizontal = _horizontalScoreOfCell(previousHorizontal) +
//                                         scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal);
//    activeCell._horizontalScore = _max(tmpGapOpenHorizontal, tmpGapExtendHorizontal);
//    activeCell._score = activeCell._horizontalScore;
//
//    if (!IsTracebackEnabled_<TTracebackConfig>::VALUE)
//        return TraceBitMap_::NONE;
//
//    TTraceValue traceGapOpen = TraceBitMap_::NONE;
//    _conditionalOrOnInequality(traceGapOpen, _horizontalScoreOfCell(activeCell), tmpGapExtendHorizontal, TraceBitMap_::HORIZONTAL_OPEN);
//    return traceGapOpen |
//           _retrieveTraceAffine(activeCell._score, TScoreValue(), tmpGapExtendHorizontal, tmpGapOpenHorizontal,
//                                TScoreValue(), TScoreValue(), RecursionDirectionHorizontal());
}

// ----------------------------------------------------------------------------
// Function _doComputeScore                        [RecursionVerticalDirection]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig>
inline typename TraceBitMap_::TTraceValue
_doComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                DPCell_<TScoreValue, AffineGaps> const & /*previousDiagonal*/,
                DPCell_<TScoreValue, AffineGaps> const & /*previousHorizontal*/,
                DPCell_<TScoreValue, AffineGaps> const & previousVertical,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionVertical const &,
                DPProfile_<TAlgorithm, AffineGaps, TTracebackConfig> const &)
{
    //typedef typename TraceBitMap_::TTraceValue TTraceValue;
    TScoreValue tmpGapOpenVertical = _scoreOfCell(previousVertical) + scoreGapOpenVertical(scoringScheme, seqHVal, seqVVal);
    activeCell._verticalScore = _verticalScoreOfCell(previousVertical) + scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal);

    // Here we distinguish between vertical and vertical open.
    activeCell._horizontalScore = DPCellDefaultInfinity<DPCell_<TScoreValue, AffineGaps> >::VALUE;
    return _internalComputeScore(activeCell, tmpGapOpenVertical, TraceBitMap_::VERTICAL, TraceBitMap_::VERTICAL_OPEN, TTracebackConfig(), RecursionDirectionVertical()) | TraceBitMap_::MAX_FROM_VERTICAL_MATRIX;
//    typedef typename TraceBitMap_::TTraceValue TTraceValue;
//    TScoreValue tmpGapOpenVertical = _scoreOfCell(previousVertical) +
//                                     scoreGapOpenVertical(scoringScheme, seqHVal, seqVVal);
//    TScoreValue tmpGapExtendVertical = _verticalScoreOfCell(previousVertical) +
//                                       scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal);
//    activeCell._verticalScore = _max(tmpGapExtendVertical, tmpGapOpenVertical);
//    activeCell._score = activeCell._verticalScore;
//
//    if (!IsTracebackEnabled_<TTracebackConfig>::VALUE)
//        return TraceBitMap_::NONE;
//
//    TTraceValue traceGapOpen = TraceBitMap_::NONE;
//    _conditionalOrOnInequality(traceGapOpen, _verticalScoreOfCell(activeCell), tmpGapExtendVertical, TraceBitMap_::VERTICAL_OPEN);
//    return traceGapOpen |
//           _retrieveTraceAffine(activeCell._score, TScoreValue(), TScoreValue(), TScoreValue(), tmpGapExtendVertical,
//                                tmpGapOpenVertical, RecursionDirectionVertical());
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_FORMULA_AFFINE_H_
