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
// Defines the recursion formula for the dp-alignment algorithms.
// ==========================================================================

// TODO(holtgrew): Documentation in this header necessary or internal only?

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag RecursionDirectionDiagonal
// ----------------------------------------------------------------------------

struct RecursionDirectionDiagonal_;
typedef Tag<RecursionDirectionDiagonal_> RecursionDirectionDiagonal;

// ----------------------------------------------------------------------------
// Tag RecursionDirectionHorizontal
// ----------------------------------------------------------------------------

struct RecursionDirectionHorizontal_;
typedef Tag<RecursionDirectionHorizontal_> RecursionDirectionHorizontal;

// ----------------------------------------------------------------------------
// Tag RecursionDirectionVertical
// ----------------------------------------------------------------------------

struct RecursionDirectionVertical_;
typedef Tag<RecursionDirectionVertical_> RecursionDirectionVertical;

// ----------------------------------------------------------------------------
// Tag RecursionDirectionAll
// ----------------------------------------------------------------------------

struct RecursionDirectionAll_;
typedef Tag<RecursionDirectionAll_> RecursionDirectionAll;

// ----------------------------------------------------------------------------
// Tag RecursionDirectionUpperDiagonal
// ----------------------------------------------------------------------------

struct RecursionDirectionUpperDiagonal_;
typedef Tag<RecursionDirectionUpperDiagonal_> RecursionDirectionUpperDiagonal;

// ----------------------------------------------------------------------------
// Tag RecursionDirectionLowerDiagonal
// ----------------------------------------------------------------------------

struct RecursionDirectionLowerDiagonal_;
typedef Tag<RecursionDirectionLowerDiagonal_> RecursionDirectionLowerDiagonal;

// ----------------------------------------------------------------------------
// Tag RecursionDirectionZero
// ----------------------------------------------------------------------------

struct RecursionDirectionZero_;
typedef Tag<RecursionDirectionZero_> RecursionDirectionZero;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _conditionalOrOnEquality()
// ----------------------------------------------------------------------------

// Function used to compare two trace values and to add a given state to the result
// value if they are equal using a bit-or operation.
template <typename TTraceValue, typename TScoreValue>
inline void
_conditionalOrOnEquality(TTraceValue & target,
                         TScoreValue const & leftComp,
                         TScoreValue const & rightComp,
                         TTraceValue state)
{
    if (leftComp == rightComp)
        target |= state;
}

// ----------------------------------------------------------------------------
// Function _conditionalOrOnInequality()
// ----------------------------------------------------------------------------

// Function used to compare two trace values and to add a given state to the result
// value if they are equal using a bit-or operation.
template <typename TTraceValue, typename TScoreValue>
inline void
_conditionalOrOnInequality(TTraceValue & target,
                           TScoreValue const & leftComp,
                           TScoreValue const & rightComp,
                           TTraceValue state)
{
    if (leftComp != rightComp)
        target |= state;
}

// ----------------------------------------------------------------------------
// Function _computeScore
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TGapCosts, typename TSequenceHValue, typename TSequenceVValue,
          typename TScoringScheme, typename TRecursionDirection, typename TDPProfile>
inline typename TraceBitMap_::TTraceValue
_computeScore(DPCell_<TScoreValue, TGapCosts> & activeCell,
              DPCell_<TScoreValue, TGapCosts> const & previousDiagonal,
              DPCell_<TScoreValue, TGapCosts> const & previousHorizontal,
              DPCell_<TScoreValue, TGapCosts> const & previousVertical,
              TSequenceHValue const & seqHVal,
              TSequenceVValue const & seqVVal,
              TScoringScheme const & scoringScheme,
              TRecursionDirection const & recDir,
              TDPProfile const & dpProfile)
{
    typedef typename TraceBitMap_::TTraceValue TTraceValue;

    TTraceValue traceDir = _doComputeScore(activeCell, previousDiagonal, previousHorizontal, previousVertical, seqHVal,
                                           seqVVal, scoringScheme, recDir, dpProfile);
    if (IsLocalAlignment_<TDPProfile>::VALUE)
        if (activeCell._score <= 0)
        {
            _setScoreOfCell(activeCell, static_cast<TScoreValue>(0));
            _setHorizontalScoreOfCell(activeCell, static_cast<TScoreValue>(0));
            _setVerticalScoreOfCell(activeCell, static_cast<TScoreValue>(0));
            return TraceBitMap_::NONE;
        }

    return traceDir;
}

// ----------------------------------------------------------------------------
// Function _doComputeScore                        [RecursionDirectionDiagonal]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TGapCosts, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TDPProfile>
inline typename TraceBitMap_::TTraceValue
_doComputeScore(DPCell_<TScoreValue, TGapCosts> & activeCell,
                DPCell_<TScoreValue, TGapCosts> const & previousDiagonal,
                DPCell_<TScoreValue, TGapCosts> const & /*previousHorizontal*/,
                DPCell_<TScoreValue, TGapCosts> const & /*previousVertical*/,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionDiagonal const &,
                TDPProfile const &)
{
    activeCell._score = _scoreOfCell(previousDiagonal) + score(scoringScheme, seqHVal, seqVVal);
    setGapExtension(activeCell, False(), False());
    if (!IsTracebackEnabled_<TDPProfile>::VALUE)
        return TraceBitMap_::NONE;

    return TraceBitMap_::DIAGONAL;
}

// ----------------------------------------------------------------------------
// Function _doComputeScore                            [RecursionDirectionZero]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TGapCosts, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgoTag, typename TTraceFlag>
inline typename TraceBitMap_::TTraceValue
_doComputeScore(DPCell_<TScoreValue, TGapCosts> & activeCell,
                DPCell_<TScoreValue, TGapCosts> const & /*previousDiagonal*/,
                DPCell_<TScoreValue, TGapCosts> const & /*previousHorizontal*/,
                DPCell_<TScoreValue, TGapCosts> const & /*previousVertical*/,
                TSequenceHValue const & /*seqHVal*/,
                TSequenceVValue const & /*seqVVal*/,
                TScoringScheme const & /*scoringScheme*/,
                RecursionDirectionZero const &,
                DPProfile_<TAlgoTag, TGapCosts, TTraceFlag> const &)
{
    _scoreOfCell(activeCell) = 0;
    return TraceBitMap_::NONE;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_H_
