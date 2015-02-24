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
// Implements the tracback algorithm for the banded chain alignment.
// First the section where both dp matrices intersect is computed to fing
// the glue poit between the traces of the current matrix and the next one.
// Afterwards the final traceback beginnging in this glue point is computed
// fir the current dp matrix.
// ==========================================================================

#ifndef INCLUDE_SEQAN_SEEDS_BANDED_CHAIN_ALIGNMENT_TRACEBACK_H_
#define INCLUDE_SEQAN_SEEDS_BANDED_CHAIN_ALIGNMENT_TRACEBACK_H_

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

// ----------------------------------------------------------------------------
// Metafunction PreferGapsAtEnd_
// ----------------------------------------------------------------------------

template <typename TFreeEndGaps, typename TMatrixSpec, typename TTracebackSpec>
struct PreferGapsAtEnd_<DPProfile_<BandedChainAlignment_<TFreeEndGaps, TMatrixSpec>,
                                   AffineGaps, TTracebackSpec> > : False{};

template <typename TFreeEndGaps, typename TTracebackSpec>
struct PreferGapsAtEnd_<DPProfile_<BandedChainAlignment_<TFreeEndGaps, BandedChainFinalDPMatrix>,
                                   AffineGaps, TTracebackSpec> > : True{};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _adaptLocalTracesToGlobalGrid()
// ----------------------------------------------------------------------------

template <typename TTraceSet, typename TGridPoint>
inline void
_adaptLocalTracesToGlobalGrid(TTraceSet & traceSet, TGridPoint const & gridBegin)
{
    typedef typename Value<TTraceSet>::Type TTraceString;
    typedef typename Iterator<TTraceString>::Type TTraceStringIterator;

    for (unsigned i = 0; i < length(traceSet); ++i)
    {
        for (TTraceStringIterator it = begin(traceSet[i]); it != end(traceSet[i]); ++it)
        {
            value(it)._horizontalBeginPos += gridBegin.i1;
            value(it)._verticalBeginPos += gridBegin.i2;
        }
    }
}

// ----------------------------------------------------------------------------
// Function _smoothGluePoint()
// ----------------------------------------------------------------------------

template <typename TTraceSegment, typename TStringSpec, typename TSize>
inline void
_smoothGluePoint(String<TTraceSegment, TStringSpec> & tracePath, TSize referenceSize)
{
    typedef typename Iterator<String<TTraceSegment, TStringSpec> >::Type TTraceSegmentsIterator;

    TTraceSegmentsIterator itEndOldTrace = end(tracePath) - referenceSize;
    TTraceSegmentsIterator itBeginNewTrace = itEndOldTrace - 1;
    if (_getTraceValue(value(itEndOldTrace)) == _getTraceValue(value(itBeginNewTrace)))
    {
        _setLength(value(itEndOldTrace), length(value(itEndOldTrace)) + length(value(itBeginNewTrace)));
        erase(tracePath, position(itBeginNewTrace, tracePath));
    }
}

// ----------------------------------------------------------------------------
// Function _glueTracebacks()
// ----------------------------------------------------------------------------

// TODO(rmaerker): What about using a set to find matching parts more efficiently. Optimize this code here!
template <typename TTraceSet>
inline void _glueTracebacks(TTraceSet & globalTraces, TTraceSet & localTraces)
{
    typedef typename Value<TTraceSet>::Type TTraceSegments;
    typedef typename Value<TTraceSegments>::Type TTraceSegment;
    typedef typename Size<TTraceSegments>::Type TSize;

    bool isGlued = false;

    TSize lengthGlobalTraces = length(globalTraces);
    TSize oldNumOfGlobalTraces = lengthGlobalTraces;
    String<unsigned> elementsToErase;

    for (unsigned j = 0; j < lengthGlobalTraces; ++j)
    {
        // traceback from back to front -> first trace segment is at end of sequences
        TTraceSegment globalTraceEndPoint = value(begin(globalTraces[j]));
        TSize numOfCurrElements = length(globalTraces[j]);
        TSize numOfAddedTraces = 0;
        // check for all existing paths if traces can be glued together.
        bool isConnected = false;
        for (unsigned i = 0; i < length(localTraces); ++i)
        {
            // traceback from back to front -> last trace segment is at beginning of sequences.
            TTraceSegment localTraceBeginPoint = value(end(localTraces[i]) - 1);
            // trace segments match in horizontal position
            if (_getEndHorizontal(globalTraceEndPoint) == _getBeginHorizontal(localTraceBeginPoint))
            {
                // trace segments match in vertical position
                if (_getEndVertical(globalTraceEndPoint) == _getBeginVertical(localTraceBeginPoint))
                { // found a glue point between local and global trace
                    TSize numOfElementsToAdd = length(localTraces[i]);

                    if (isConnected)
                    {
                        // create a new traceback track.
                        TTraceSegments newTraceTrack;
                        resize(newTraceTrack, numOfCurrElements + numOfElementsToAdd, Generous());
                        arrayMoveForward(begin(localTraces[i]), end(localTraces[i]), begin(newTraceTrack));
                        arrayCopyForward(end(globalTraces[j]) - numOfCurrElements, end(globalTraces[j]), begin(newTraceTrack) + numOfElementsToAdd);
                        appendValue(globalTraces, newTraceTrack);
                        ++numOfAddedTraces;
                        continue;
                    }
                    // resize global traces such that new elements fit into it
                    resize(globalTraces[j], numOfCurrElements + numOfElementsToAdd, Generous());
                    // shift old values to correct position in array after new elements will be added.
                    // use backward move in order to avoid overwriting when ranges overlap
                    arrayMoveBackward(begin(globalTraces[j]), begin(globalTraces[j]) + numOfCurrElements,
                                      begin(globalTraces[j]) + numOfElementsToAdd);
                    // actually move new elements to global traces at it's begin
                    arrayMoveForward(begin(localTraces[i]), end(localTraces[i]), begin(globalTraces[j]));
                    isGlued = true;
                    isConnected = true;
                }
            }
        }
        if (!isConnected)
        {
            appendValue(elementsToErase, j);
        }
        else
        {
            // First, smooth the actual path we are currently gluing.
            _smoothGluePoint(globalTraces[j], numOfCurrElements);
            // Second, smooth all paths that have been added to the global traceback.
            for (unsigned traceId = oldNumOfGlobalTraces; traceId < oldNumOfGlobalTraces + numOfAddedTraces; ++traceId)
                _smoothGluePoint(globalTraces[traceId], numOfCurrElements);
            oldNumOfGlobalTraces += numOfAddedTraces;
        }
    }

    for (unsigned i = length(elementsToErase); i > 0; --i)
    {
        erase(globalTraces, elementsToErase[i-1]); // erase from behind to avoid accessing an element beyond the scope
    }
    SEQAN_ASSERT_EQ_MSG(isGlued, true, "Fatal error while trying to connect trace backs: No glue point available!");
}

// ----------------------------------------------------------------------------
// Function _correctDPCellForAffineGaps()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TTraceValue>
inline void
_correctDPCellForAffineGaps(DPCell_<TScoreValue, LinearGaps> const &, TTraceValue /*lastTraceValue*/)
{
    //no-op
}

template <typename TScoreValue, typename TTraceValue>
inline void
_correctDPCellForAffineGaps(DPCell_<TScoreValue, AffineGaps> & dpCell, TTraceValue lastTraceValue)
{
    typedef DPCell_<TScoreValue, AffineGaps> TDPCell;
    if (lastTraceValue & TraceBitMap_::DIAGONAL)
    {
        dpCell._verticalScore = DPCellDefaultInfinity<TDPCell>::VALUE;
        dpCell._horizontalScore = DPCellDefaultInfinity<TDPCell>::VALUE;
    }
    else if (lastTraceValue & TraceBitMap_::VERTICAL)
        dpCell._horizontalScore = DPCellDefaultInfinity<TDPCell>::VALUE;
    else
        dpCell._verticalScore = DPCellDefaultInfinity<TDPCell>::VALUE;
}

// ----------------------------------------------------------------------------
// Function _computeTraceback()                          [BandedChainAlignment]
// ----------------------------------------------------------------------------

template<typename TTarget, typename TDPTraceMatrixNavigator, typename TDPCell, typename TScoutSpec,
         typename TSequenceH, typename TSequenceV, typename TBandFlag, typename TFreeEndGaps, typename TDPMatrixLocation,
         typename TGapCosts, typename TTracebackSpec>
void _computeTraceback(TTarget & target,
                       TDPTraceMatrixNavigator & matrixNavigator,
                       unsigned maxHostPosition,
                       DPScout_<TDPCell, TScoutSpec> & dpScout,
                       TSequenceH const & seqH,
                       TSequenceV const & seqV,
                       DPBandConfig<TBandFlag> const & band,
                       DPProfile_<BandedChainAlignment_<TFreeEndGaps, TDPMatrixLocation>, TGapCosts, TTracebackSpec> const & dpProfile)
{
    typedef DPScout_<TDPCell, TScoutSpec> TDPScout_;
    typedef typename TDPScout_::TScoutState TScoutState_;
    typedef typename TScoutState_::TInitCell TInitCell;
    typedef typename Container<TDPTraceMatrixNavigator>::Type TContainer;
    typedef typename Size<TContainer>::Type TSize;
    //typedef typename MakeSigned<TSize>::Type TSignedSize;
    typedef typename Position<TContainer>::Type TPosition;
    typedef typename MakeSigned<TPosition>::Type TSignedPosition;
    typedef typename Size<TSequenceH>::Type TSizeSeqH;
    typedef typename Size<TSequenceV>::Type TSizeSeqV;
    typedef typename TraceBitMap_::TTraceValue TTraceValue;

    if (IsSameType<TTracebackSpec, TracebackOff>::VALUE)
        return;

    TSizeSeqH seqHSize = length(seqH);
    TSizeSeqV seqVSize = length(seqV);

    // Determine whether or not we place gaps to the left.
    typedef typename IsGapsLeft_<TTracebackSpec>::Type TIsGapsLeft;

    // TODO(rmaerker): Define separate function for this.
    // Set to the correct position within the trace matrix.
    _setToPosition(matrixNavigator, maxHostPosition);

    SEQAN_ASSERT_LEQ(coordinate(matrixNavigator, +DPMatrixDimension_::HORIZONTAL), seqHSize);
    SEQAN_ASSERT_LEQ(coordinate(matrixNavigator, +DPMatrixDimension_::VERTICAL), seqVSize);

    TTraceValue traceValue = value(matrixNavigator);
    TTraceValue lastTraceValue = _retrieveInitialTraceDirection(traceValue, dpProfile);

    TracebackCoordinator_<TPosition> tracebackCoordinator(coordinate(matrixNavigator, +DPMatrixDimension_::HORIZONTAL),
                                                          coordinate(matrixNavigator, +DPMatrixDimension_::VERTICAL),
                                                          dpScout._dpScoutStatePtr->_horizontalNextGridOrigin,
                                                          dpScout._dpScoutStatePtr->_verticalNextGridOrigin,
                                                          band, seqHSize, seqVSize);

    // Record trailing gaps if any.
    if (IsSameType<TDPMatrixLocation, BandedChainFinalDPMatrix>::VALUE)
    {
        if (tracebackCoordinator._currRow != seqVSize)
            _recordSegment(target, seqHSize, tracebackCoordinator._currRow, seqVSize - tracebackCoordinator._currRow,
                           +TraceBitMap_::VERTICAL);
        if (tracebackCoordinator._currColumn != seqHSize)
            _recordSegment(target, tracebackCoordinator._currColumn, tracebackCoordinator._currRow, seqHSize -
                           tracebackCoordinator._currColumn, +TraceBitMap_::HORIZONTAL);

        _computeTraceback(target, matrixNavigator, position(matrixNavigator), seqH, seqV, band, dpProfile);
        return;
    }

    TSize fragmentLength = 0;
    TTarget tmp;
    while(!_hasReachedEnd(tracebackCoordinator) && traceValue != +TraceBitMap_::NONE)
        _doTraceback(tmp, matrixNavigator, traceValue, lastTraceValue, fragmentLength, tracebackCoordinator, TGapCosts(), TIsGapsLeft());

    TSignedPosition horizontalInitPos = static_cast<TSignedPosition>(tracebackCoordinator._currColumn) -
                                        static_cast<TSignedPosition>(tracebackCoordinator._endColumn);
    TSignedPosition verticalInitPos = static_cast<TSignedPosition>(tracebackCoordinator._currRow) -
                                      static_cast<TSignedPosition>(tracebackCoordinator._endRow);

    bool insertResult;
    if (verticalInitPos <= 0)
    {
        _correctDPCellForAffineGaps(dpScout._dpScoutStatePtr->_horizontalInitNextMatrix[horizontalInitPos], lastTraceValue);
        insertResult = dpScout._dpScoutStatePtr->_nextInitializationCells.insert(TInitCell(horizontalInitPos, 0,
                                                dpScout._dpScoutStatePtr->_horizontalInitNextMatrix[horizontalInitPos])).second;
    }
    else
    {
        _correctDPCellForAffineGaps(dpScout._dpScoutStatePtr->_verticalInitNextMatrix[verticalInitPos], lastTraceValue);
        insertResult = dpScout._dpScoutStatePtr->_nextInitializationCells.insert(TInitCell(0, verticalInitPos,
                                                dpScout._dpScoutStatePtr->_verticalInitNextMatrix[verticalInitPos])).second;
    }

    // Now before we can continue at the current position, we might want to track a vertical/horizontal gap up to this position.
    if (insertResult)
    {
        if (verticalInitPos < 0)  // Here we are in a vertical gap.
            _recordSegment(target, tracebackCoordinator._currColumn, tracebackCoordinator._currRow, -verticalInitPos,
                           lastTraceValue);
        else if (horizontalInitPos < 0)  // Here we are in a horizontal gap.
            _recordSegment(target, tracebackCoordinator._currColumn, tracebackCoordinator._currRow, -horizontalInitPos,
                           lastTraceValue);
        _computeTraceback(target, matrixNavigator, position(matrixNavigator), seqH, seqV, band, dpProfile);
    }

    if (IsSameType<TDPMatrixLocation, BandedChainInitialDPMatrix>::VALUE)
    {
        TPosition currCol = coordinate(matrixNavigator, +DPMatrixDimension_::HORIZONTAL);
        TPosition currRow = coordinate(matrixNavigator, +DPMatrixDimension_::VERTICAL);

        // Correct the row position.
        if (IsSameType<TBandFlag, BandOn>::VALUE)
            if (upperDiagonal(band) > 0)
                if (currCol < tracebackCoordinator._breakpoint1)
                    if (currCol < tracebackCoordinator._breakpoint2)
                        currRow -= length(container(matrixNavigator), +DPMatrixDimension_::VERTICAL) - 1 + lowerDiagonal(band) - currCol;
        // Record leading gaps if any.
        if (currRow != 0u)
            _recordSegment(target, 0, 0, currRow, +TraceBitMap_::VERTICAL);
        if (currCol != 0u)
            _recordSegment(target, 0, 0, currCol, +TraceBitMap_::HORIZONTAL);
    }
}

// ----------------------------------------------------------------------------
// Function _computeTraceback()        [BandedChainAlignment, global interface]
// ----------------------------------------------------------------------------

template <typename TTarget, typename TDPTraceMatrixNavigator, typename TDPScout, typename TSequenceH, typename TSequenceV,
typename TBandFlag, typename TFreeEndGaps, typename TDPMatrixLocation, typename TGapCosts, typename TTracebackSpec>
void _computeTraceback(StringSet<TTarget> & targetSet,
                       TDPTraceMatrixNavigator & matrixNavigator,
                       TDPScout & dpScout,
                       TSequenceH const & seqH,
                       TSequenceV const & seqV,
                       DPBandConfig<TBandFlag> const & band,
                       DPProfile_<BandedChainAlignment_<TFreeEndGaps, TDPMatrixLocation>, TGapCosts, TTracebackSpec> const & dpProfile)
{
    typedef typename TDPScout::TMaxHostPositionString TMaxHostPositions;
    typedef typename Iterator<TMaxHostPositions, Standard>::Type TMaxHostPositionsIterator;

    // We have to clear the nextInitalization cells here.
    dpScout._dpScoutStatePtr->_nextInitializationCells.clear();
    TMaxHostPositions & tracebackCandidates = maxHostPositions(dpScout);

    TMaxHostPositionsIterator itTraceCandidates = begin(tracebackCandidates, Standard());

    for (;itTraceCandidates != end(tracebackCandidates, Standard()); ++itTraceCandidates)
    {
        TTarget tmpTarget;
        _computeTraceback(tmpTarget, matrixNavigator, *itTraceCandidates, dpScout, seqH, seqV, band, dpProfile);
        if (!empty(tmpTarget))
            appendValue(targetSet, tmpTarget);  // We only need to add the tmpTarget if it gives us a new alignment.
    }
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_SEEDS_BANDED_CHAIN_ALIGNMENT_TRACEBACK_H_
