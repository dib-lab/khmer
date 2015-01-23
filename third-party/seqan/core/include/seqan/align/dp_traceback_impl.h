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
// Implements the traceback algorithm.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_TRACEBACK_IMPL_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_TRACEBACK_IMPL_H_

// TODO(holtgrew): GapsRight traceback is currently untested.
// TODO(rmaerker): Change Tracback to TraceConfig<TGapsPlacement, TPathSpec> | TraceBackOff

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class TracebackCoordinator_
// ----------------------------------------------------------------------------

template <typename TPosition>
class TracebackCoordinator_
{
public:
    TPosition _currColumn;
    TPosition _currRow;
    TPosition _endColumn;
    TPosition _endRow;
    TPosition _breakpoint1;      // First breakpoint where banded trace switches to unbanded trace.
    TPosition _breakpoint2;      // Second breakpoint where unbanded trace switches back to banded trace. Only if begin of upper diagonal is bigger than end of lower diagonal.
    bool _isInBand;

    template <typename TBandFlag, typename TSizeH, typename TSizeV>
    TracebackCoordinator_(TPosition currColumn,
                          TPosition currRow,
                          DPBand_<TBandFlag> const & band,
                          TSizeH seqHSize,
                          TSizeV seqVSize)
        : _currColumn(currColumn),
          _currRow(currRow),
          _endColumn(0u),
          _endRow(0u),
          _breakpoint1(0u),
          _breakpoint2(0u),
          _isInBand(false)
    {
        _initTracebackCoordinator(*this, band, seqHSize, seqVSize);
    }

    template <typename TBandFlag, typename TSizeH, typename TSizeV>
    TracebackCoordinator_(TPosition currColumn,
                          TPosition currRow,
                          TPosition endColumn,
                          TPosition endRow,
                          DPBand_<TBandFlag> const & band,
                          TSizeH seqHSize,
                          TSizeV seqVSize)
        : _currColumn(currColumn),
          _currRow(currRow),
          _endColumn(endColumn),
          _endRow(endRow),
          _breakpoint1(0u),
          _breakpoint2(0u),
          _isInBand(false)
    {
        _initTracebackCoordinator(*this, band, seqHSize, seqVSize);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction PreferGapsAtEnd_
// ----------------------------------------------------------------------------

// Checks whether the gaps at the end should be preferred over a matching area.
template <typename TDPProfile>
struct PreferGapsAtEnd_ : False{};

template <typename TAlgorithm, typename TTracebackSpec>
struct PreferGapsAtEnd_<DPProfile_<TAlgorithm, AffineGaps, TTracebackSpec > > : True{};

template <typename TAlgorithm, typename TTraceSpec>
struct PreferGapsAtEnd_<DPProfile_<TAlgorithm, LinearGaps, TracebackOn<TracebackConfig_<TTraceSpec, GapsRight> > > > : True{};


// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _hasReachedEnd()
// ----------------------------------------------------------------------------

template <typename TPosition>
inline bool
_hasReachedEnd(TracebackCoordinator_<TPosition> const & coordinator)
{
    return coordinator._currColumn <= coordinator._endColumn || coordinator._currRow <= coordinator._endRow;
}

// ----------------------------------------------------------------------------
// Function _initTracebackCoordinator()
// ----------------------------------------------------------------------------

template <typename TPosition, typename TBandFlag, typename TSizeH, typename TSizeV>
inline void
_initTracebackCoordinator(TracebackCoordinator_<TPosition> & coordinator,
                          DPBand_<TBandFlag> const & band,
                          TSizeH seqHSize,
                          TSizeV seqVSize)
{
    typedef typename Position<DPBand_<TBandFlag> >::Type TBandPosition;
    if (IsSameType<TBandFlag, BandOn>::VALUE)
    {
        // Adapt the current column value when the lower diagonal is positive (shift right in horizontal direction).
        if (lowerDiagonal(band) >= 0)
            coordinator._currColumn += static_cast<TPosition>(lowerDiagonal(band));
        // Adapt the current row value when the current column comes after the upper diagonal (shift down in vertical direction).
        if (static_cast<TBandPosition>(coordinator._currColumn) > upperDiagonal(band))
            coordinator._currRow += coordinator._currColumn - upperDiagonal(band);
        // Adapt the end row value when the end column comes after the upper diagonal (shift down in vertical direction).
        if (static_cast<TBandPosition>(coordinator._endColumn) > upperDiagonal(band))
            coordinator._endRow += coordinator._endColumn - upperDiagonal(band);

        coordinator._breakpoint1 = _min(seqHSize, static_cast<TSizeH>(_max(0, upperDiagonal(band))));
        coordinator._breakpoint2 = _min(seqHSize, static_cast<TSizeH>(_max(0, static_cast<TBandPosition>(seqVSize) +
                                                                           lowerDiagonal(band))));
        // Update the current row if the current column is before the upper diagoal or the first column where the maximal band size is reached.
        if (coordinator._currColumn < _min(coordinator._breakpoint1, coordinator._breakpoint2))
            coordinator._currRow -= _min(coordinator._breakpoint1, coordinator._breakpoint2) - coordinator._currColumn;
        coordinator._isInBand = true;
    }
}

// ----------------------------------------------------------------------------
// Function _isInBand()
// ----------------------------------------------------------------------------

template <typename TPosition>
inline bool
_isInBand(TracebackCoordinator_<TPosition> const & coordinator)
{
    if (!coordinator._isInBand)
        return coordinator._isInBand;
    return (coordinator._currColumn > coordinator._breakpoint1 || coordinator._currColumn <= coordinator._breakpoint2);
}


// ----------------------------------------------------------------------------
// Function _doTracebackGoDiagonal()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TDPTraceMatrixNavigator, typename TTraceValue, typename TSize, typename TPosition,
          typename TGapCosts>
inline void
_doTracebackGoDiagonal(TTarget & target,
                       TDPTraceMatrixNavigator & matrixNavigator,
                       TTraceValue & traceValue,
                       TTraceValue & lastTraceValue,
                       TSize & fragmentLength,
                       TracebackCoordinator_<TPosition> & tracebackCoordinator,
                       TGapCosts const &)
{
    if (!(lastTraceValue & TraceBitMap_::DIAGONAL)) // the old trace value was not diagonal
    {
        _recordSegment(target, tracebackCoordinator._currColumn, tracebackCoordinator._currRow, fragmentLength,
                       lastTraceValue);
        
        lastTraceValue = TraceBitMap_::DIAGONAL;
        fragmentLength = 0;
    }
    _traceDiagonal(matrixNavigator, _isInBand(tracebackCoordinator));
    traceValue = value(matrixNavigator);
    --tracebackCoordinator._currColumn;
    --tracebackCoordinator._currRow;
    ++fragmentLength;
}

// ----------------------------------------------------------------------------
// Function _doTracebackGoVertical()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TDPTraceMatrixNavigator, typename TTraceValue, typename TSize, typename TPosition,
          typename TGapCosts>
inline void
_doTracebackGoVertical(TTarget & target,
                       TDPTraceMatrixNavigator & matrixNavigator,
                       TTraceValue & traceValue,
                       TTraceValue & lastTraceValue,
                       TSize & fragmentLength,
                       TracebackCoordinator_<TPosition> & tracebackCoordinator,
                       TGapCosts const &)
{
    if (!(lastTraceValue & TraceBitMap_::VERTICAL)) // the old trace value was not diagonal
    {
        _recordSegment(target, tracebackCoordinator._currColumn, tracebackCoordinator._currRow, fragmentLength,
                       lastTraceValue);

        lastTraceValue = TraceBitMap_::VERTICAL;
        fragmentLength = 0;
    }
    // We are in a vertical gap. So continue after we reach the end of the vertical gap.
    if (IsSameType<TGapCosts, AffineGaps>::VALUE)
    {
        while ((!(traceValue & TraceBitMap_::VERTICAL_OPEN) || (traceValue & TraceBitMap_::VERTICAL)) && (tracebackCoordinator._currRow != 1))
        {
            _traceVertical(matrixNavigator, _isInBand(tracebackCoordinator));
            traceValue = value(matrixNavigator);
            --tracebackCoordinator._currRow;
            ++fragmentLength;
        }
        // We have to ensure, that we do not continue in vertical direction if we reached a vertical_open sign.
        _traceVertical(matrixNavigator, _isInBand(tracebackCoordinator));
        // Forbid continuing in vertical direction.
        traceValue = value(matrixNavigator); // & (TraceBitMap_::NO_VERTICAL_TRACEBACK | TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX);
        --tracebackCoordinator._currRow;
        ++fragmentLength;
    }
    else
    {
        _traceVertical(matrixNavigator, _isInBand(tracebackCoordinator));
        traceValue = value(matrixNavigator);
        --tracebackCoordinator._currRow;
        ++fragmentLength;
    }
}

// ----------------------------------------------------------------------------
// Function _doTracebackMaxFromVertical()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TDPTraceMatrixNavigator, typename TTraceValue, typename TSize, typename TPosition,
          typename TGapCosts>
inline void
_doTracebackMaxFromVertical(TTarget & target,
                            TDPTraceMatrixNavigator & matrixNavigator,
                            TTraceValue & traceValue,
                            TTraceValue & lastTraceValue,
                            TSize & fragmentLength,
                            TracebackCoordinator_<TPosition> & tracebackCoordinator,
                            TGapCosts const &)
{
    if (!(lastTraceValue & TraceBitMap_::VERTICAL)) // the old trace value was not diagonal
    {
        _recordSegment(target, tracebackCoordinator._currColumn, tracebackCoordinator._currRow, fragmentLength,
                       lastTraceValue);
        lastTraceValue = TraceBitMap_::VERTICAL;
        fragmentLength = 0;
    }
    _traceVertical(matrixNavigator, _isInBand(tracebackCoordinator));
    // Forbid continuing in vertical direction.
    traceValue = value(matrixNavigator); // & (TraceBitMap_::NO_VERTICAL_TRACEBACK | TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX);
    --tracebackCoordinator._currRow;
    ++fragmentLength;
}

// ----------------------------------------------------------------------------
// Function _doTracebackGoHorizontal()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TDPTraceMatrixNavigator, typename TTraceValue, typename TSize, typename TPosition,
          typename TGapCosts>
inline void
_doTracebackGoHorizontal(TTarget & target,
                         TDPTraceMatrixNavigator & matrixNavigator,
                         TTraceValue & traceValue,
                         TTraceValue & lastTraceValue,
                         TSize & fragmentLength,
                         TracebackCoordinator_<TPosition> & tracebackCoordinator,
                         TGapCosts const &)
{
    if (!(lastTraceValue & TraceBitMap_::HORIZONTAL)) // the old trace value was not diagonal
    {
        _recordSegment(target, tracebackCoordinator._currColumn, tracebackCoordinator._currRow, fragmentLength,
                       lastTraceValue);

        lastTraceValue = TraceBitMap_::HORIZONTAL;
        fragmentLength = 0;
    }
    if (IsSameType<TGapCosts, AffineGaps>::VALUE)
    {
        while ((!(traceValue & TraceBitMap_::HORIZONTAL_OPEN) || (traceValue & TraceBitMap_::HORIZONTAL)) && (tracebackCoordinator._currColumn != 1))
        {
            _traceHorizontal(matrixNavigator, _isInBand(tracebackCoordinator));
            traceValue = value(matrixNavigator);
            --tracebackCoordinator._currColumn;
            ++fragmentLength;
        }
        _traceHorizontal(matrixNavigator, _isInBand(tracebackCoordinator));
        // Forbid continuing in horizontal direction.
        traceValue = value(matrixNavigator);  // & (TraceBitMap_::NO_HORIZONTAL_TRACEBACK | TraceBitMap_::MAX_FROM_VERTICAL_MATRIX);
        --tracebackCoordinator._currColumn;
        ++fragmentLength;
    }
    else
    {
        _traceHorizontal(matrixNavigator, _isInBand(tracebackCoordinator));
        traceValue = value(matrixNavigator);
        --tracebackCoordinator._currColumn;
        ++fragmentLength;
    }
}

// ----------------------------------------------------------------------------
// Function _doTracebackMaxFromHorizontal()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TDPTraceMatrixNavigator, typename TTraceValue, typename TSize, typename TPosition,
          typename TGapCosts>
inline void
_doTracebackMaxFromHorizontal(TTarget & target,
                              TDPTraceMatrixNavigator & matrixNavigator,
                              TTraceValue & traceValue,
                              TTraceValue & lastTraceValue,
                              TSize & fragmentLength,
                              TracebackCoordinator_<TPosition> & tracebackCoordinator,
                              TGapCosts const &)
{
    if (!(lastTraceValue & TraceBitMap_::HORIZONTAL)) // the old trace value was not diagonal
    {
        _recordSegment(target, tracebackCoordinator._currColumn, tracebackCoordinator._currRow, fragmentLength,
                       lastTraceValue);
        lastTraceValue = TraceBitMap_::HORIZONTAL;
        fragmentLength = 0;
    }
    _traceHorizontal(matrixNavigator, _isInBand(tracebackCoordinator));
    // Forbid continuing in horizontal direction.
    traceValue = value(matrixNavigator); // & (TraceBitMap_::NO_HORIZONTAL_TRACEBACK | TraceBitMap_::MAX_FROM_VERTICAL_MATRIX);
    --tracebackCoordinator._currColumn;
    ++fragmentLength;
}

// ----------------------------------------------------------------------------
// Function _doTraceback()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TDPTraceMatrixNavigator, typename TTraceValue, typename TSize, typename TPosition,
          typename TGapCosts, typename TIsGapsLeft>
inline void
_doTraceback(TTarget & target,
             TDPTraceMatrixNavigator & matrixNavigator,
             TTraceValue & traceValue,
             TTraceValue & lastTraceValue,
             TSize & fragmentLength,
             TracebackCoordinator_<TPosition> & tracebackCoordinator,
             TGapCosts const & gapsCost,
             TIsGapsLeft const & /*isGapsLeft*/)
{
    if (TIsGapsLeft::VALUE)  // Gaps should be placed on the left.
    {
        if (traceValue & TraceBitMap_::DIAGONAL)
        {
            _doTracebackGoDiagonal(target, matrixNavigator, traceValue, lastTraceValue, fragmentLength, tracebackCoordinator, gapsCost);
        }  // In case of Gotoh we prefer the longest possible way in this direction.
        else if (traceValue & TraceBitMap_::MAX_FROM_VERTICAL_MATRIX && traceValue & TraceBitMap_::VERTICAL)
        {
            _doTracebackGoVertical(target, matrixNavigator, traceValue, lastTraceValue, fragmentLength, tracebackCoordinator, gapsCost);
        }
        else if (traceValue & TraceBitMap_::MAX_FROM_VERTICAL_MATRIX && traceValue & TraceBitMap_::VERTICAL_OPEN)
        {
            _doTracebackMaxFromVertical(target, matrixNavigator, traceValue, lastTraceValue, fragmentLength, tracebackCoordinator, gapsCost);
        }
        else if (traceValue & TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX && traceValue & TraceBitMap_::HORIZONTAL)
        {
            _doTracebackGoHorizontal(target, matrixNavigator, traceValue, lastTraceValue, fragmentLength, tracebackCoordinator, gapsCost);
        }
        else if (traceValue & TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX && traceValue & TraceBitMap_::HORIZONTAL_OPEN)
        {
            _doTracebackMaxFromHorizontal(target, matrixNavigator, traceValue, lastTraceValue, fragmentLength, tracebackCoordinator, gapsCost);
        }  // In case of Gotoh we prefer the longest possible way in this direction.
        else // the trace back is either NONE or something else
        {
            if (traceValue == TraceBitMap_::NONE)
            {
                return;
            }
            SEQAN_ASSERT_FAIL("Reached undefined traceback value!");
        }
    }
    else  // Gaps should be placed on the right.
    {
        if (traceValue & TraceBitMap_::MAX_FROM_VERTICAL_MATRIX && traceValue & TraceBitMap_::VERTICAL)
        {
            _doTracebackGoVertical(target, matrixNavigator, traceValue, lastTraceValue, fragmentLength, tracebackCoordinator, gapsCost);
        }
        else if (traceValue & TraceBitMap_::MAX_FROM_VERTICAL_MATRIX && traceValue & TraceBitMap_::VERTICAL_OPEN)
        {
            _doTracebackMaxFromVertical(target, matrixNavigator, traceValue, lastTraceValue, fragmentLength, tracebackCoordinator, gapsCost);
        }
        else if (traceValue & TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX && traceValue & TraceBitMap_::HORIZONTAL)
        {
            _doTracebackGoHorizontal(target, matrixNavigator, traceValue, lastTraceValue, fragmentLength, tracebackCoordinator, gapsCost);
        }
        else if (traceValue & TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX && traceValue & TraceBitMap_::HORIZONTAL_OPEN)
        {
            _doTracebackMaxFromHorizontal(target, matrixNavigator, traceValue, lastTraceValue, fragmentLength, tracebackCoordinator, gapsCost);
        }  // In case of Gotoh we prefer the longest possible way in this direction.
        else if (traceValue & TraceBitMap_::DIAGONAL)
        {
            _doTracebackGoDiagonal(target, matrixNavigator, traceValue, lastTraceValue, fragmentLength, tracebackCoordinator, gapsCost);
        }  // In case of Gotoh we prefer the longest possible way in this direction.
        else // the trace back is either NONE or something else
        {
            if (traceValue == TraceBitMap_::NONE)
            {
                return;
            }
            SEQAN_ASSERT_FAIL("Reached undefined traceback value!");
        }
    }
}

// ----------------------------------------------------------------------------
// Function _retrieveInitialTraceDirection()
// ----------------------------------------------------------------------------

template <typename TTraceValue, typename TDPProfile>
inline TTraceValue
_retrieveInitialTraceDirection(TTraceValue & traceValue, TDPProfile const & /*dpProfile*/)
{
    if (PreferGapsAtEnd_<TDPProfile>::VALUE)
    {
        if (traceValue & TraceBitMap_::MAX_FROM_VERTICAL_MATRIX)
        {
            traceValue &= (TraceBitMap_::VERTICAL | TraceBitMap_::VERTICAL_OPEN | TraceBitMap_::MAX_FROM_VERTICAL_MATRIX);
            return TraceBitMap_::VERTICAL;
        }
        else if (traceValue & TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX)
        {
            traceValue &= (TraceBitMap_::HORIZONTAL | TraceBitMap_::HORIZONTAL_OPEN | TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX);
            return TraceBitMap_::HORIZONTAL;
        }
        return TraceBitMap_::DIAGONAL;  // We set the last value to the 
    }

    if (traceValue & TraceBitMap_::DIAGONAL)
            return TraceBitMap_::DIAGONAL;
    if (traceValue & (TraceBitMap_::VERTICAL | TraceBitMap_::MAX_FROM_VERTICAL_MATRIX))
        return  TraceBitMap_::VERTICAL;
    if (traceValue & (TraceBitMap_::HORIZONTAL | TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX))
        return TraceBitMap_::HORIZONTAL;

    return TraceBitMap_::NONE;
}

// ----------------------------------------------------------------------------
// Function _computeTraceback()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TDPTraceMatrixNavigator, typename TSequenceH, typename TSequenceV,
          typename TBandFlag, typename TAlgorithm, typename TGapCosts, typename TTracebackSpec>
void _computeTraceback(TTarget & target,
                       TDPTraceMatrixNavigator & matrixNavigator,
                       unsigned  maxHostPosition,
                       TSequenceH const & seqH,
                       TSequenceV const & seqV,
                       DPBand_<TBandFlag> const & band,
                       DPProfile_<TAlgorithm, TGapCosts, TTracebackSpec> const & dpProfile)
{
    typedef typename Container<TDPTraceMatrixNavigator>::Type TContainer;
    typedef typename Size<TContainer>::Type TSize;
    typedef typename Position<TContainer>::Type TPosition;
    typedef typename TraceBitMap_::TTraceValue TTraceValue;
    typedef typename Size<TSequenceH>::Type TSizeH;
    typedef typename Size<TSequenceV>::Type TSizeV;

    if (IsSameType<TTracebackSpec, TracebackOff>::VALUE)
        return;

    // Determine whether or not we place gaps to the left.
    typedef typename IsGapsLeft_<TTracebackSpec>::Type TIsGapsLeft;

    TSizeH seqHSize = length(seqH);
    TSizeV seqVSize = length(seqV);

    // Set the navigator to the position where the maximum was found.
    _setToPosition(matrixNavigator, maxHostPosition);

    SEQAN_ASSERT_LEQ(coordinate(matrixNavigator, +DPMatrixDimension_::HORIZONTAL), seqHSize);
    SEQAN_ASSERT_LEQ(coordinate(matrixNavigator, +DPMatrixDimension_::VERTICAL), seqVSize);

    TTraceValue traceValue = value(matrixNavigator);
    TTraceValue lastTraceValue = _retrieveInitialTraceDirection(traceValue, dpProfile);

    TracebackCoordinator_<TPosition> tracebackCoordinator(coordinate(matrixNavigator, +DPMatrixDimension_::HORIZONTAL),
                                                          coordinate(matrixNavigator, +DPMatrixDimension_::VERTICAL),
                                                          band, seqHSize, seqVSize);

    if (IsGlobalAlignment_<TAlgorithm>::VALUE)
    {
        if (tracebackCoordinator._currRow != seqVSize)
            _recordSegment(target, seqHSize, tracebackCoordinator._currRow, seqVSize - tracebackCoordinator._currRow,
                           +TraceBitMap_::VERTICAL);
        if (tracebackCoordinator._currColumn != seqHSize)
            _recordSegment(target, tracebackCoordinator._currColumn, tracebackCoordinator._currRow, seqHSize -
                           tracebackCoordinator._currColumn, +TraceBitMap_::HORIZONTAL);
    }

    TSize fragmentLength = 0;
    while (!_hasReachedEnd(tracebackCoordinator) && traceValue != TraceBitMap_::NONE)
        _doTraceback(target, matrixNavigator, traceValue, lastTraceValue, fragmentLength, tracebackCoordinator, TGapCosts(), TIsGapsLeft());


    // Record last detected fragment.
    _recordSegment(target, tracebackCoordinator._currColumn, tracebackCoordinator._currRow, fragmentLength, lastTraceValue);
    if (IsGlobalAlignment_<TAlgorithm>::VALUE)
    {
        // Record leading gaps if any.
        if (tracebackCoordinator._currRow != 0u)
            _recordSegment(target, 0, 0, tracebackCoordinator._currRow, +TraceBitMap_::VERTICAL);
        if (tracebackCoordinator._currColumn != 0u)
            _recordSegment(target, 0, 0, tracebackCoordinator._currColumn, +TraceBitMap_::HORIZONTAL);
    }
}

// Needed as a delegation method to allow invocation of both methods with host position and dpScout.
template <typename TTarget, typename TDPTraceMatrixNavigator, typename TDPCell, typename TScoutSpec,
          typename TSequenceH, typename TSequenceV, typename TBandFlag, typename TAlgorithm, typename TGapCosts,
          typename TTracebackSpec>
void _computeTraceback(TTarget & target,
                       TDPTraceMatrixNavigator & matrixNavigator,
                       DPScout_<TDPCell, TScoutSpec> const & dpScout,
                       TSequenceH const & seqH,
                       TSequenceV const & seqV,
                       DPBand_<TBandFlag> const & band,
                       DPProfile_<TAlgorithm, TGapCosts, TTracebackSpec> const & dpProfile)
{
    _computeTraceback(target, matrixNavigator, maxHostPosition(dpScout), seqH, seqV, band, dpProfile);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_TRACEBACK_IMPL_H_
