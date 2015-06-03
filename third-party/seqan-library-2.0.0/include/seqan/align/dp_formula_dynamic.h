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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Implements the score function for dynamic gap costs published in
// "Dynamic Gaps Selector: A Smith Waterman Sequence Alignment Algorithm with
// Affine Gap Model Optimization" by Gianvito Urgese et al.
// ==========================================================================

#ifndef INCLUDE_SEQAN_ALIGN_DP_FORMULA_DYNAMIC_H_
#define INCLUDE_SEQAN_ALIGN_DP_FORMULA_DYNAMIC_H_

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
// Function _internalComputeScore      [RecursionDirectionDiagonal, DynamicGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueGap>
inline typename TraceBitMap_::TTraceValue
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      TScoreValue const & diagCompare,
                      TTraceValueL,
                      TTraceValueGap,
                      TracebackOff const &,
                      RecursionDirectionDiagonal const &)
{
    if(_scoreOfCell(activeCell) < diagCompare)
    {
        activeCell._score = diagCompare;
        setGapExtension(activeCell, False(), False());
        return TraceBitMap_::NONE;
    }
    return TraceBitMap_::NONE;
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueGap, typename TGapsPlacement>
inline typename TraceBitMap_::TTraceValue
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      TScoreValue const & diagCompare,
                      TTraceValueL leftTrace,
                      TTraceValueGap gapTrace,
                      TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> > const &,
                      RecursionDirectionDiagonal const &)
{
    if(_scoreOfCell(activeCell) <= diagCompare)
    {
        activeCell._score = diagCompare;
        setGapExtension(activeCell, False(), False());
        return TraceBitMap_::DIAGONAL | leftTrace;
    }
    return leftTrace | gapTrace;
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueGap, typename TGapsPlacement>
inline typename TraceBitMap_::TTraceValue
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      TScoreValue const & diagCompare,
                      TTraceValueL leftTrace,
                      TTraceValueGap gapTrace,
                      TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> >  const &,
                      RecursionDirectionDiagonal const &)
{
    if(_scoreOfCell(activeCell) < diagCompare)
    {
        activeCell._score = diagCompare;  // Maximal score comes from diagonal.
        setGapExtension(activeCell, False(), False());
        return TraceBitMap_::DIAGONAL | leftTrace;  // Return trace for Diagonal.
    }
    if (_scoreOfCell(activeCell) == diagCompare)  // Maximal score comes from previous computed directions and diagonal.
        return leftTrace | TraceBitMap_::DIAGONAL | gapTrace;  // Return all directions inclusively the flag indicating max from gap.

    return leftTrace | gapTrace;  // Maximum comes from gap. Return gap value inclusively the flag indicating max from gap.
}

// ----------------------------------------------------------------------------
// Function _internalComputeScore    [RecursionDirectionHorizontal, DynamicGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TValueH, typename TValueV, typename TScore>
inline typename TraceBitMap_::TTraceValue
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      DPCell_<TScoreValue, DynamicGaps> const & prevCell,
                      TValueH const & valH,
                      TValueV const & valV,
                      TScore const & score,
                      TracebackOff const &,
                      RecursionDirectionHorizontal const &)
{
    if(!isGapExtension(prevCell, DynamicGapExtensionHorizontal()))
        activeCell._score = _scoreOfCell(prevCell) + scoreGapOpenHorizontal(score, valH, valV);
    return TraceBitMap_::NONE;
}

template <typename TScoreValue, typename TValueH, typename TValueV, typename TScore, typename TTraceConfig>
inline typename TraceBitMap_::TTraceValue
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      DPCell_<TScoreValue, DynamicGaps> const & prevCell,
                      TValueH const & valH,
                      TValueV const & valV,
                      TScore const & score,
                      TracebackOn<TTraceConfig> const &,
                      RecursionDirectionHorizontal const &)
{
    if (!isGapExtension(prevCell, DynamicGapExtensionHorizontal()))
    {
        activeCell._score = _scoreOfCell(prevCell) + scoreGapOpenHorizontal(score, valH, valV);
        return TraceBitMap_::HORIZONTAL_OPEN;
    }
    return TraceBitMap_::HORIZONTAL;
}

// ----------------------------------------------------------------------------
// Function _internalComputeScore      [RecursionDirectionVertical, DynamicGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TValueH, typename TValueV, typename TScore>
inline typename TraceBitMap_::TTraceValue
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      DPCell_<TScoreValue, DynamicGaps> const & prevCell,
                      TValueH const & valH,
                      TValueV const & valV,
                      TScore const & score,
                      TracebackOff const &,
                      RecursionDirectionVertical const &)
{
    if(!isGapExtension(prevCell, DynamicGapExtensionVertical()))
        activeCell._score = _scoreOfCell(prevCell) + scoreGapOpenVertical(score, valH, valV);
    return TraceBitMap_::NONE;
}

template <typename TScoreValue, typename TValueH, typename TValueV, typename TScore, typename TTraceConfig>
inline typename TraceBitMap_::TTraceValue
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      DPCell_<TScoreValue, DynamicGaps> const & prevCell,
                      TValueH const & valH,
                      TValueV const & valV,
                      TScore const & score,
                      TracebackOn<TTraceConfig> const &,
                      RecursionDirectionVertical const &)
{
    if (!isGapExtension(prevCell, DynamicGapExtensionVertical()))
    {
        activeCell._score = _scoreOfCell(prevCell) + scoreGapOpenVertical(score, valH, valV);
        return TraceBitMap_::VERTICAL_OPEN;
    }
    return TraceBitMap_::VERTICAL;
}

// ----------------------------------------------------------------------------
// Function _internalComputeScore          [Vertical vs Horizontal, DynamicGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue>
inline typename TraceBitMap_::TTraceValue
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      TScoreValue const & horizontalComp,
                      TracebackOff const &)
{
    if(_scoreOfCell(activeCell) < horizontalComp)
    {
        activeCell._score = horizontalComp;
        setGapExtension(activeCell, False(), True());
        return TraceBitMap_::NONE;
    }
    setGapExtension(activeCell, True(), False());
    return TraceBitMap_::NONE;
}

template <typename TScoreValue, typename TGapsPlacement>
inline typename TraceBitMap_::TTraceValue
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      TScoreValue const & horizontalComp,
                      TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> >  const &)
{
    if(_scoreOfCell(activeCell) < horizontalComp)
    {
        activeCell._score = horizontalComp;
        setGapExtension(activeCell, False(), True());
        return TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX;
    }
    setGapExtension(activeCell, True(), False());
    return TraceBitMap_::MAX_FROM_VERTICAL_MATRIX;
}

template <typename TScoreValue, typename TGapsPlacement>
inline typename TraceBitMap_::TTraceValue
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      TScoreValue const & horizontalComp,
                      TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> >  const &)
{
    if(_scoreOfCell(activeCell) < horizontalComp)
    {
        setGapExtension(activeCell, False(), True());
        activeCell._score = horizontalComp;
        return TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX;
    }
    if (_scoreOfCell(activeCell) == horizontalComp)
    {
        setGapExtension(activeCell, True(), True());
        return TraceBitMap_::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX;
    }
    setGapExtension(activeCell, True(), False());
    return TraceBitMap_::MAX_FROM_VERTICAL_MATRIX;
}

// ----------------------------------------------------------------------------
// Function _doComputeScore                 [RecursionAllDirection, DynamicGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig>
inline typename TraceBitMap_::TTraceValue
_doComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                DPCell_<TScoreValue, DynamicGaps> const & previousDiagonal,
                DPCell_<TScoreValue, DynamicGaps> const & previousHorizontal,
                DPCell_<TScoreValue, DynamicGaps> const & previousVertical,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionAll const &,
                DPProfile_<TAlgorithm, DynamicGaps, TTracebackConfig> const &)
{
    typedef typename TraceBitMap_::TTraceValue TTraceValue;

    // Compute best alignment from either horizontal open or extension.
    DPCell_<TScoreValue, DynamicGaps> tmpScore = _scoreOfCell(previousHorizontal) + scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal);
    TTraceValue tvGap = _internalComputeScore(tmpScore, previousHorizontal, seqHVal, seqVVal, scoringScheme,
                                              TTracebackConfig(), RecursionDirectionHorizontal());

    // Compute best alignment between vertical and vertical open gap.
    activeCell._score = _scoreOfCell(previousVertical) + scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal);
    tvGap |= _internalComputeScore(activeCell, previousVertical, seqHVal, seqVVal, scoringScheme,
                                   TTracebackConfig(), RecursionDirectionVertical());

    // Finds the maximum between the vertical and the horizontal matrix. Stores the flag for coming from a potential direction.
    TTraceValue tvMax = _internalComputeScore(activeCell, tmpScore._score, TTracebackConfig());  // Stores from where the maximal score comes.
    tmpScore._score = _scoreOfCell(previousDiagonal) + score(scoringScheme, seqHVal, seqVVal);
    return _internalComputeScore(activeCell, tmpScore._score, tvGap, tvMax, TTracebackConfig(), RecursionDirectionDiagonal());
}

// ----------------------------------------------------------------------------
// Function _doComputeScore       [RecursionUpperDiagonalDirection, DynamicGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig>
inline typename TraceBitMap_::TTraceValue
_doComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                DPCell_<TScoreValue, DynamicGaps> const & previousDiagonal,
                DPCell_<TScoreValue, DynamicGaps> const & previousHorizontal,
                DPCell_<TScoreValue, DynamicGaps> const & /*previousVertical*/,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionUpperDiagonal const &,
                DPProfile_<TAlgorithm, DynamicGaps, TTracebackConfig> const &)
{
    typedef typename TraceBitMap_::TTraceValue TTraceValue;

    // This computes the difference between the horizontal extend and horizontal open.
    activeCell._score = _scoreOfCell(previousHorizontal) + scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal);
    TTraceValue tv = _internalComputeScore(activeCell, previousHorizontal, seqHVal, seqVVal, scoringScheme,
                                           TTracebackConfig(), RecursionDirectionHorizontal());

    setGapExtension(activeCell, False(), True());
    TScoreValue tmpScore = _scoreOfCell(previousDiagonal) + score(scoringScheme, seqHVal, seqVVal);
    return _internalComputeScore(activeCell, tmpScore, tv, TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX,
                                 TTracebackConfig(), RecursionDirectionDiagonal());
}

// ----------------------------------------------------------------------------
// Function _doComputeScore      [RecursionDirectionLowerDiagonal, DynamicGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig>
inline typename TraceBitMap_::TTraceValue
_doComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                DPCell_<TScoreValue, DynamicGaps> const & previousDiagonal,
                DPCell_<TScoreValue, DynamicGaps> const & /*previousHorizontal*/,
                DPCell_<TScoreValue, DynamicGaps> const & previousVertical,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionLowerDiagonal const &,
                DPProfile_<TAlgorithm, DynamicGaps, TTracebackConfig> const &)
{
    typedef typename TraceBitMap_::TTraceValue TTraceValue;

    // This computes the difference between the vertical extend and vertical open.
    activeCell._score = _scoreOfCell(previousVertical) + scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal);
    TTraceValue tv = _internalComputeScore(activeCell, previousVertical, seqHVal, seqVVal, scoringScheme,
                                           TTracebackConfig(), RecursionDirectionVertical());
    setGapExtension(activeCell, True(), False());
    TScoreValue tmpScore = _scoreOfCell(previousDiagonal) + score(scoringScheme, seqHVal, seqVVal);
    return _internalComputeScore(activeCell, tmpScore, tv, TraceBitMap_::MAX_FROM_VERTICAL_MATRIX,
                                 TTracebackConfig(), RecursionDirectionDiagonal());
}

// ----------------------------------------------------------------------------
// Function _doComputeScore                      [RecursionHorizontalDirection]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig>
inline typename TraceBitMap_::TTraceValue
_doComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                DPCell_<TScoreValue, DynamicGaps> const & /*previousDiagonal*/,
                DPCell_<TScoreValue, DynamicGaps> const & previousHorizontal,
                DPCell_<TScoreValue, DynamicGaps> const & /*previousVertical*/,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionHorizontal const & tag,
                DPProfile_<TAlgorithm, DynamicGaps, TTracebackConfig> const &)
{
    activeCell._score = _scoreOfCell(previousHorizontal) + scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal);
    setGapExtension(activeCell, False(), True());
    return _internalComputeScore(activeCell, previousHorizontal, seqHVal, seqVVal, scoringScheme,
                                 TTracebackConfig(), tag) | TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX;
}

// ----------------------------------------------------------------------------
// Function _doComputeScore                        [RecursionVerticalDirection]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig>
inline typename TraceBitMap_::TTraceValue
_doComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                DPCell_<TScoreValue, DynamicGaps> const & /*previousDiagonal*/,
                DPCell_<TScoreValue, DynamicGaps> const & /*previousHorizontal*/,
                DPCell_<TScoreValue, DynamicGaps> const & previousVertical,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionVertical const & tag,
                DPProfile_<TAlgorithm, DynamicGaps, TTracebackConfig> const &)
{
    activeCell._score = _scoreOfCell(previousVertical) + scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal);
    setGapExtension(activeCell, True(), False());
    return _internalComputeScore(activeCell, previousVertical, seqHVal, seqVVal, scoringScheme,
                                 TTracebackConfig(), tag) | TraceBitMap_::MAX_FROM_VERTICAL_MATRIX;
}

}  // namespace seqan

#endif  // INCLUDE_SEQAN_ALIGN_DP_FORMULA_DYNAMIC_H_
