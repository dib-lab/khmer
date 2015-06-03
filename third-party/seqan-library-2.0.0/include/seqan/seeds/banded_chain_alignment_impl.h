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
// Implementation of a global alignment between two sequences given a
// set of seeds in monotonic non-decreasing order. This implementation is
// based on the algorithm described in the "LAGAN and Multi-LAGAN: Efficient
// Tools for Large-Scale Multiple Alignment of Genomic Data", Brudno et al.,
// Genome Res. 2003, 13: 721-731, doi:10.1101/gr.926603
// ==========================================================================

#ifndef INCLUDE_SEQAN_SEEDS_BANDED_CHAIN_ALIGNMENT_IMPL_H_
#define INCLUDE_SEQAN_SEEDS_BANDED_CHAIN_ALIGNMENT_IMPL_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Struct BandedChainTracking()
// ----------------------------------------------------------------------------

// Used to distinguish between a cell that has to be tracked for the
// initialization of the next grid and the cell that has to be tracked to find
// the optimal score for the current traceback.

struct BandedChainTracking
{
    enum Options
    {
        OPTION_INIT = 0,
        OPTION_IS_LAST_COLUMN = 1,
        OPTION_IS_LAST_ROW = 2,
        OPTION_STORE_INIT_COLUMN = 4,
        OPTION_STORE_INIT_ROW = 8
    };
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _isLastSeed()
// ----------------------------------------------------------------------------

template <typename TSeedSet>
inline bool
_isLastSeed(typename Iterator<TSeedSet const, Standard>::Type iter,
            TSeedSet const & seedSet)
{
    typedef typename Iterator<TSeedSet const>::Type TSeedIterator;
    TSeedIterator itEnd = end(seedSet, Standard());
    return iter == --itEnd;
}

// ----------------------------------------------------------------------------
// Function _checkColinearity()
// ----------------------------------------------------------------------------

template <typename TIter>
inline bool
_checkColinearity(TIter currIter)
{
    TIter nextIter = currIter;
    ++nextIter;

    return endPositionH(value(currIter)) <= beginPositionH(value(nextIter)) &&
           endPositionV(value(currIter)) <= beginPositionV(value(nextIter));
}

// ----------------------------------------------------------------------------
// Function _isCrossingBeginHorizontal()
// ----------------------------------------------------------------------------

// Checks whether a given seed crosses the first row of the global grid.
template <typename TPosition, typename TSize>
inline bool _isCrossingBeginHorizontal(TPosition const & horizontalSeedBeginPos, TSize const & bandExtension)
{
    typedef typename MakeSigned<TSize>::Type TSignedSize;
    if ((TSignedSize) horizontalSeedBeginPos - (TSignedSize) bandExtension <= 0)
        return true;
    return false;
}

// ----------------------------------------------------------------------------
// Function _isCrossingBeginVertical()
// ----------------------------------------------------------------------------

// Checks whether a given seed crosses the first column of the global grid.
template <typename TPosition, typename TSize>
inline bool _isCrossingBeginVertical(TPosition const & verticalSeedBeginPos, TSize const & bandExtension)
{
    typedef typename MakeSigned<TSize>::Type TSignedSize;

    if ((TSignedSize) verticalSeedBeginPos - (TSignedSize) bandExtension <= 0)
        return true;
    return false;
}

// ----------------------------------------------------------------------------
// Function _isCrossingEndHorizontal()
// ----------------------------------------------------------------------------

// Checks whether a given seed crosses the last row of the global grid.
template <typename TPosition, typename TSize, typename TSequenceH>
inline bool
_isCrossingEndHorizontal(TPosition const & horizontalSeedEndPos,
                         TSize const & bandSize,
                         TSequenceH const & seqH)
{
    typedef typename MakeSigned<TSize>::Type TSignedSize;

    if ((TSignedSize) horizontalSeedEndPos + (TSignedSize) bandSize >= (TSignedSize) length(seqH))
        return true;
    return false;
}

// ----------------------------------------------------------------------------
// Function _isCrossingEndVertical()
// ----------------------------------------------------------------------------

// Checks whether a given seed crosses the last column of the global grid.
template <typename TPosition, typename TSize, typename TSequenceV>
inline bool
_isCrossingEndVertical(TPosition const & verticalSeedEndPos,
                       TSize const & bandSize,
                       TSequenceV const & seqV)
{
    typedef typename MakeSigned<TSize>::Type TSignedSize;

    if ((TSignedSize) verticalSeedEndPos + (TSignedSize) bandSize >= (TSignedSize) length(seqV))
        return true;
    return false;
}

// ----------------------------------------------------------------------------
// Function _horizontalBandShiftBeginPoint()
// ----------------------------------------------------------------------------

// Returns the shift in horizontal direction for skew bands.
template <typename TSeed>
inline typename Size<TSeed>::Type _horizontalBandShiftBeginPoint(TSeed const & seed)
{
    typedef typename Size<TSeed>::Type TSize;
    TSize bandSize = (upperDiagonal(seed) - (beginPositionH(seed) - beginPositionV(seed)));
    SEQAN_ASSERT_GEQ(bandSize, TSize(0));  // must be greater or equal than zero
    return bandSize;
}

// ----------------------------------------------------------------------------
// Function _verticalBandShiftBeginPoint()
// ----------------------------------------------------------------------------

// Returns the shift in vertical direction for skew bands.
template <typename TSeed>
inline typename Size<TSeed>::Type _verticalBandShiftBeginPoint(TSeed const & seed)
{
    typedef typename Size<TSeed>::Type TSize;
    TSize bandSize((beginPositionH(seed) - beginPositionV(seed)) - lowerDiagonal(seed));
    SEQAN_ASSERT_GEQ(bandSize, TSize(0));  // must be greater or equal than zero
    return bandSize;
}

// ----------------------------------------------------------------------------
// Function _horizontalBandShiftEndPoint()
// ----------------------------------------------------------------------------

    // Returns the shift in horizontal direction for skew bands.
template <typename TSeed>
inline typename Size<TSeed>::Type _horizontalBandShiftEndPoint(TSeed const & seed)
{
    typedef typename Size<TSeed>::Type TSize;
    TSize bandSize = endPositionH(seed) - endPositionV(seed) - lowerDiagonal(seed);
    SEQAN_ASSERT_GEQ(bandSize, TSize(0));  // must be greater or equal than zero
    return bandSize;
}

// ----------------------------------------------------------------------------
// Function _verticalBandShiftEndPoint()
// ----------------------------------------------------------------------------

    // Returns the shift in vertical direction for skew bands.
template <typename TSeed>
inline typename Size<TSeed>::Type _verticalBandShiftEndPoint(TSeed const & seed)
{
    typedef typename Size<TSeed>::Type TSize;
    TSize bandSize = upperDiagonal(seed) - endPositionH(seed) + endPositionV(seed);
    SEQAN_ASSERT_GEQ(bandSize, TSize(0));  // must be greater or equal than zero
    return bandSize;
}


// ----------------------------------------------------------------------------
// Function _verticalCellInitialization()
// ----------------------------------------------------------------------------

// Initiliazes vertical cell with corresponding value from dp scout state.
template <typename TDPCell, typename TTraceNavigator>
inline TDPCell &
_verticalCellInitialization(DPScout_<TDPCell, BandedChainAlignmentScout> const & dpScout,
                            TTraceNavigator const & navigator)
{
    return dpScout._dpScoutStatePtr->_verticalInitCurrentMatrix[coordinate(navigator, +DPMatrixDimension_::VERTICAL) -
                                                                 (navigator._laneLeap - 1)];
}

// ----------------------------------------------------------------------------
// Function _horizontalCellInitialization()
// ----------------------------------------------------------------------------

// Initializes horizontal cell with corresponding value from dp scout state.
template <typename TDPCell, typename TTraceNavigator>
inline TDPCell const &
_horizontalCellInitialization(DPScout_<TDPCell, BandedChainAlignmentScout> const & dpScout,
                              TTraceNavigator const & navigator)
{
    return dpScout._dpScoutStatePtr->_horizontalInitCurrentMatrix[coordinate(navigator, +DPMatrixDimension_::HORIZONTAL)];
}

// ----------------------------------------------------------------------------
// Function _determineTrackingOptions()
// ----------------------------------------------------------------------------

// Determines the various tracking options for the current dp matrix.
template <typename TScout, typename TNavigator, typename TColumnDescriptor, typename TCellDescriptor, typename TFreeEndGaps,
          typename TMatrixSpec, typename TGapCosts, typename TTracebackSpec>
inline void _determineTrackingOptions(unsigned & res,
                                      TScout const & scout,
                                      TNavigator const & traceMatrixNavigator,
                                      TColumnDescriptor const & /*columnDescriptor*/,
                                      TCellDescriptor const & /*cellDescriptor*/,
                                      DPProfile_<BandedChainAlignment_<TFreeEndGaps, TMatrixSpec>, TGapCosts, TTracebackSpec> const &)
{
    if (coordinate(traceMatrixNavigator, +DPMatrixDimension_::HORIZONTAL) >=
        scout._dpScoutStatePtr->_horizontalNextGridOrigin)
    {
        // Matches an initialization row.
        if (IsSameType<typename TColumnDescriptor::TLocation, PartialColumnBottom>::VALUE)
        {
            // Matches the cells that have to be stored for initializing the next matrix horizontally.
            if (coordinate(traceMatrixNavigator, +DPMatrixDimension_::VERTICAL) + traceMatrixNavigator._laneLeap ==
                scout._dpScoutStatePtr->_verticalNextGridOrigin)
                res |= BandedChainTracking::OPTION_STORE_INIT_ROW; //BANDED_CHAIN_TRACKING_OPTION_STORE_INIT_ROW;
        }  // Matches an initialization row.
        else
        {
            // Matches the cells that have to be stored for initializing the next matrix horizontally.
            if (coordinate(traceMatrixNavigator, +DPMatrixDimension_::VERTICAL) ==
                scout._dpScoutStatePtr->_verticalNextGridOrigin)
                res |= BandedChainTracking::OPTION_STORE_INIT_ROW;
        }

        // Matches an initialization column.
        if (coordinate(traceMatrixNavigator, +DPMatrixDimension_::HORIZONTAL) ==
                    scout._dpScoutStatePtr->_horizontalNextGridOrigin)
            // Matches the cells that have to be stored for initializing the next matrix vertically.
            if (coordinate(traceMatrixNavigator, +DPMatrixDimension_::VERTICAL) >=
                scout._dpScoutStatePtr->_verticalNextGridOrigin)
                res |= BandedChainTracking::OPTION_STORE_INIT_COLUMN;

        // We are in the column that has to be tracked, if we are in the last cell.
        if (IsSameType<TCellDescriptor, LastCell>::VALUE)
        {
            if (IsSameType<TMatrixSpec, BandedChainFinalDPMatrix>::VALUE)
            {
                if (IsFreeEndGap_<TFreeEndGaps, DPLastRow>::VALUE)
                    res |= BandedChainTracking::OPTION_IS_LAST_ROW;
            }
            else
            {
                res |= BandedChainTracking::OPTION_IS_LAST_ROW;
            }
        }

        // We track the maximal score if we are in the final column
        // For the full column we need to check if we are beyond the initialization row.
        if (IsSameType<typename TColumnDescriptor::TLocation, FullColumn>::VALUE)
        {
            if (IsSameType<typename TColumnDescriptor::TColumnProperty, DPFinalColumn>::VALUE)
            {
                if (IsSameType<TCellDescriptor, LastCell>::VALUE)
                {
                    res |= BandedChainTracking::OPTION_IS_LAST_COLUMN | BandedChainTracking::OPTION_IS_LAST_ROW;
                }
                else if (coordinate(traceMatrixNavigator, +DPMatrixDimension_::VERTICAL) >=
                         scout._dpScoutStatePtr->_verticalNextGridOrigin)
                {
                    if (IsSameType<TMatrixSpec, BandedChainFinalDPMatrix>::VALUE)
                    {
                        if (IsFreeEndGap_<TFreeEndGaps, DPLastColumn>::VALUE)
                            res |= BandedChainTracking::OPTION_IS_LAST_COLUMN;
                    }
                    else
                    {
                        res |= BandedChainTracking::OPTION_IS_LAST_COLUMN;
                    }
                }
            }
        }
        else  // Banded version we track all scores - initialization row ends in first cell of final column
        {
            if (IsSameType<typename TColumnDescriptor::TColumnProperty, DPFinalColumn>::VALUE)
            {
                if (IsSameType<TCellDescriptor, LastCell>::VALUE)
                {
                    res |= BandedChainTracking::OPTION_IS_LAST_COLUMN | BandedChainTracking::OPTION_IS_LAST_ROW;
                }
                else
                {
                    if (IsSameType<TMatrixSpec, BandedChainFinalDPMatrix>::VALUE)
                    {
                        if (IsFreeEndGap_<TFreeEndGaps, DPLastColumn>::VALUE)
                            res |= BandedChainTracking::OPTION_IS_LAST_COLUMN;
                    }
                    else
                    {
                        res |= BandedChainTracking::OPTION_IS_LAST_COLUMN;
                    }
                }
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function _applyBandedChainTracking()
// ----------------------------------------------------------------------------

template <typename TDPScout, typename TTraceMatrixNavigator, typename TDPCell, typename TColumnDescriptor,
          typename TCellDescriptor, typename TDPProfile>
inline void
_applyBandedChainTracking(TDPScout & scout,
                          TTraceMatrixNavigator & traceMatrixNavigator,
                          TDPCell const & activeCell,
                          TColumnDescriptor const &,
                          TCellDescriptor const &,
                          TDPProfile const &)
{
    unsigned result = BandedChainTracking::OPTION_INIT;
    _determineTrackingOptions(result, scout, traceMatrixNavigator, TColumnDescriptor(), TCellDescriptor(), TDPProfile());
    _scoutBestScore(scout, activeCell, traceMatrixNavigator,
                    (result & BandedChainTracking::OPTION_IS_LAST_COLUMN),
                    (result & BandedChainTracking::OPTION_IS_LAST_ROW),
                    (result & BandedChainTracking::OPTION_STORE_INIT_COLUMN),
                    (result & BandedChainTracking::OPTION_STORE_INIT_ROW));
}

// ----------------------------------------------------------------------------
// Function _computeCell()                               [BandedChainAlignment]
// ----------------------------------------------------------------------------

// Overload of _computeCell function to add functionality specific to the bande chain alignment.
template <typename TDPScout, typename TTraceMatrixNavigator, typename TScoreValue, typename TGapCosts,
          typename TSeqHValue, typename TSeqVValue, typename TScoringScheme, typename TColumnDescriptor,
          typename TCellDescriptor, typename TFreeEndGaps, typename TDPMatrixLocation, typename TTracebackConfig>
inline void
_computeCell(TDPScout & scout,
             TTraceMatrixNavigator & traceMatrixNavigator,
             DPCell_<TScoreValue, TGapCosts> & activeCell,
             DPCell_<TScoreValue, TGapCosts> const & previousDiagonal,
             DPCell_<TScoreValue, TGapCosts> const & previousHorizontal,
             DPCell_<TScoreValue, TGapCosts> const & previousVertical,
             TSeqHValue const & seqHVal,
             TSeqVValue const & seqVVal,
             TScoringScheme const & scoringScheme,
             TColumnDescriptor const &,
             TCellDescriptor const &,
             DPProfile_<BandedChainAlignment_<TFreeEndGaps, TDPMatrixLocation>, TGapCosts, TracebackOn<TTracebackConfig> > const &)
{
    typedef DPProfile_<BandedChainAlignment_<TFreeEndGaps, TDPMatrixLocation>, TGapCosts, TracebackOn<TTracebackConfig> > TDPProfile;
    typedef DPMetaColumn_<TDPProfile, TColumnDescriptor> TMetaColumnProfile;

    assignValue(
        traceMatrixNavigator,
        _computeScore(activeCell, previousDiagonal, previousHorizontal, previousVertical, seqHVal, seqVVal,
                      scoringScheme, typename RecursionDirection_<TMetaColumnProfile, TCellDescriptor>::Type(),
                      TDPProfile()));
    if (TrackingEnabled_<TMetaColumnProfile, TCellDescriptor>::VALUE)
        _applyBandedChainTracking(scout, traceMatrixNavigator, activeCell, TColumnDescriptor(), TCellDescriptor(), TDPProfile());
}

// ----------------------------------------------------------------------------
// Function _computeCell()              [BandedChainAlignment, DPInitialColumn]
// ----------------------------------------------------------------------------

template <typename TDPScout, typename TDPTraceMatrixNavigator, typename TScoreValue, typename TGapCosts,
          typename TSeqHValue, typename TSeqVValue, typename TScoringScheme, typename TColumnType, typename TCellDescriptor,
          typename TSpec, typename TDPMatrixLocation, typename TTracebackConfig>
inline void
_computeCell(TDPScout & scout,
             TDPTraceMatrixNavigator & traceMatrixNavigator,
             DPCell_<TScoreValue, TGapCosts> & activeCell,
             DPCell_<TScoreValue, TGapCosts> const &,
             DPCell_<TScoreValue, TGapCosts> const &,
             DPCell_<TScoreValue, TGapCosts> const &,
             TSeqHValue const &,
             TSeqVValue const &,
             TScoringScheme const & /*scoringScheme*/,
             MetaColumnDescriptor<DPInitialColumn, TColumnType> const & /*metaColumnDescriptor*/,
             TCellDescriptor const & /*cellDescriptor*/,
             DPProfile_<BandedChainAlignment_<TSpec, TDPMatrixLocation>, TGapCosts, TracebackOn<TTracebackConfig> > const & /*dpProfile*/)
{
    typedef DPProfile_<BandedChainAlignment_<TSpec, TDPMatrixLocation>, TGapCosts, TracebackOn<TTracebackConfig> > TDPProfile;
    typedef MetaColumnDescriptor<DPInitialColumn, TColumnType> TColumnDescriptor;
    typedef DPMetaColumn_<TDPProfile, TColumnDescriptor> TMetaColumnProfile;

    activeCell = _verticalCellInitialization(scout, traceMatrixNavigator);
    assignValue(traceMatrixNavigator, +TraceBitMap_::NONE);
    if (TrackingEnabled_<TMetaColumnProfile, TCellDescriptor>::VALUE)
        _applyBandedChainTracking(scout, traceMatrixNavigator, activeCell, TColumnDescriptor(), TCellDescriptor(), TDPProfile());
}

// ----------------------------------------------------------------------------
// Function _computeHorizontalInitCell()
// ----------------------------------------------------------------------------

template <typename TDPScout, typename TDPTraceMatrixNavigator, typename TDPCell,
          typename TColumnDescriptor, typename TCellDescriptor, typename TDPProfile>
inline void
_computeHorizontalInitCell(TDPScout & scout,
                           TDPTraceMatrixNavigator & traceMatrixNavigator,
                           TDPCell & activeCell,
                           TColumnDescriptor const &,
                           TCellDescriptor const &,
                           TDPProfile const &)
{
    typedef DPMetaColumn_<TDPProfile, TColumnDescriptor> TMetaColumnProfile;

    activeCell = _horizontalCellInitialization(scout, traceMatrixNavigator);
    assignValue(traceMatrixNavigator, +TraceBitMap_::NONE);
    if (TrackingEnabled_<TMetaColumnProfile, TCellDescriptor>::VALUE)
        _applyBandedChainTracking(scout, traceMatrixNavigator, activeCell, TColumnDescriptor(), TCellDescriptor(), TDPProfile());
}

// ----------------------------------------------------------------------------
// Function _computeCell()     [BandedChainAlignment, PartialColumn, FirstCell]
// ----------------------------------------------------------------------------

// For DPInnerColumn.
template <typename TDPScout, typename TDPTraceMatrixNavigator, typename TScoreValue, typename TGapCosts,
          typename TSeqHValue, typename TSeqVValue, typename TScoringScheme, typename TSpec, typename TDPMatrixLocation,
          typename TTracebackConfig>
inline void
_computeCell(TDPScout & scout,
             TDPTraceMatrixNavigator & traceMatrixNavigator,
             DPCell_<TScoreValue, TGapCosts> & activeCell,
             DPCell_<TScoreValue, TGapCosts> const &,
             DPCell_<TScoreValue, TGapCosts> const &,
             DPCell_<TScoreValue, TGapCosts> const &,
             TSeqHValue const &,
             TSeqVValue const &,
             TScoringScheme const &,
             MetaColumnDescriptor<DPInnerColumn, PartialColumnTop> const &,
             FirstCell const &,
             DPProfile_<BandedChainAlignment_<TSpec, TDPMatrixLocation>, TGapCosts, TracebackOn<TTracebackConfig> > const &)
{
    typedef DPProfile_<BandedChainAlignment_<TSpec, TDPMatrixLocation>, TGapCosts, TracebackOn<TTracebackConfig> > TDPProfile;
    _computeHorizontalInitCell(scout, traceMatrixNavigator, activeCell,
                               MetaColumnDescriptor<DPInnerColumn, PartialColumnTop>(), FirstCell(), TDPProfile());
}

// For DPFinalColumn.
template <typename TDPScout, typename TDPTraceMatrixNavigator, typename TScoreValue, typename TGapCosts,
          typename TSeqHValue, typename TSeqVValue, typename TScoringScheme, typename TSpec, typename TDPMatrixLocation,
          typename TTracebackConfig>
inline void
_computeCell(TDPScout & scout,
             TDPTraceMatrixNavigator & traceMatrixNavigator,
             DPCell_<TScoreValue, TGapCosts> & activeCell,
             DPCell_<TScoreValue, TGapCosts> const &,
             DPCell_<TScoreValue, TGapCosts> const &,
             DPCell_<TScoreValue, TGapCosts> const &,
             TSeqHValue const &,
             TSeqVValue const &,
             TScoringScheme const &,
             MetaColumnDescriptor<DPFinalColumn, PartialColumnTop> const &,
             FirstCell const &,
             DPProfile_<BandedChainAlignment_<TSpec, TDPMatrixLocation>, TGapCosts, TracebackOn<TTracebackConfig> > const &)
{
    typedef DPProfile_<BandedChainAlignment_<TSpec, TDPMatrixLocation>, TGapCosts, TracebackOn<TTracebackConfig> > TDPProfile;
    _computeHorizontalInitCell(scout, traceMatrixNavigator, activeCell,
                               MetaColumnDescriptor<DPFinalColumn, PartialColumnTop>(), FirstCell(), TDPProfile());
}

// ----------------------------------------------------------------------------
// Function _computeCell()          [BandedChainAlignment, FullColumn, FirstCell]
// ----------------------------------------------------------------------------

// For DPInnerColumn.
template <typename TDPScout, typename TDPTraceMatrixNavigator, typename TScoreValue, typename TGapCosts,
          typename TSeqHValue, typename TSeqVValue, typename TScoringScheme, typename TSpec, typename TDPMatrixLocation,
          typename TTracebackConfig>
inline void
_computeCell(TDPScout & scout,
             TDPTraceMatrixNavigator & traceMatrixNavigator,
             DPCell_<TScoreValue, TGapCosts> & activeCell,
             DPCell_<TScoreValue, TGapCosts> const &,
             DPCell_<TScoreValue, TGapCosts> const &,
             DPCell_<TScoreValue, TGapCosts> const &,
             TSeqHValue const &,
             TSeqVValue const &,
             TScoringScheme const &,
             MetaColumnDescriptor<DPInnerColumn, FullColumn> const &,
             FirstCell const &,
             DPProfile_<BandedChainAlignment_<TSpec, TDPMatrixLocation>, TGapCosts, TracebackOn<TTracebackConfig> > const &)
{
    typedef DPProfile_<BandedChainAlignment_<TSpec, TDPMatrixLocation>, TGapCosts, TracebackOn<TTracebackConfig> > TDPProfile;
    _computeHorizontalInitCell(scout, traceMatrixNavigator, activeCell,
                               MetaColumnDescriptor<DPInnerColumn, FullColumn>(), FirstCell(), TDPProfile());
}

// For DPFinalColumn.
template <typename TDPScout, typename TDPTraceMatrixNavigator, typename TScoreValue, typename TGapCosts,
          typename TSeqHValue, typename TSeqVValue, typename TScoringScheme, typename TSpec, typename TDPMatrixLocation,
          typename TTracebackConfig>
inline void
_computeCell(TDPScout & scout,
             TDPTraceMatrixNavigator & traceMatrixNavigator,
             DPCell_<TScoreValue, TGapCosts> & activeCell,
             DPCell_<TScoreValue, TGapCosts> const &,
             DPCell_<TScoreValue, TGapCosts> const &,
             DPCell_<TScoreValue, TGapCosts> const &,
             TSeqHValue const &,
             TSeqVValue const &,
             TScoringScheme const &,
             MetaColumnDescriptor<DPFinalColumn, FullColumn> const &,
             FirstCell const &,
             DPProfile_<BandedChainAlignment_<TSpec, TDPMatrixLocation>, TGapCosts, TracebackOn<TTracebackConfig> > const &)
{
    typedef DPProfile_<BandedChainAlignment_<TSpec, TDPMatrixLocation>, TGapCosts, TracebackOn<TTracebackConfig> > TDPProfile;
    _computeHorizontalInitCell(scout, traceMatrixNavigator, activeCell,
                               MetaColumnDescriptor<DPFinalColumn, FullColumn>(), FirstCell(), TDPProfile());
}

// ----------------------------------------------------------------------------
// Function findFirstAnchor()
// ----------------------------------------------------------------------------

// Finds the last seed that overlaps with the first row or column of the matrix.
// This way we guarantee that special conditions need to be considered only for one seed.
template <typename TSeedSet>
inline typename Iterator<TSeedSet const, Standard>::Type
_findFirstAnchor(TSeedSet const & seedSet, int bandExtension)
{
    typedef typename Iterator<TSeedSet const, Standard>::Type TIterator;
    typedef typename Value<TSeedSet const>::Type TSeed;

    SEQAN_ASSERT_GT_MSG(length(seedSet), 0u, "SeedSet is empty!");

    TIterator it = begin(seedSet, Standard());
    TIterator itEnd = end(seedSet, Standard());
    --itEnd;

    while (it != itEnd)
    {
        TSeed seed = value(++it);
        if (_isCrossingBeginHorizontal(beginPositionH(seed), bandExtension))
            continue;
        else if (_isCrossingBeginVertical(beginPositionV(seed), bandExtension))
            continue;
        else
        { // found seed which is not crossing the begin.
            SEQAN_ASSERT(it != static_cast<TIterator>(begin(seedSet, Standard())));
            return --it;
        }
    }
    return it;
}

// ----------------------------------------------------------------------------
// Function findLastAnchor()
// ----------------------------------------------------------------------------

// Finds the last seed that is not overlapping with the last row and column.
// This way we guarantee that special conditions need to be considered only for one seed.
template <typename TSeedSet, typename TSequenceH, typename TSequenceV>
inline typename Iterator<TSeedSet const, Standard>::Type
_findLastAnchor(typename Iterator<TSeedSet const, Standard>::Type & iterBegin,
                TSeedSet const & seedSet, TSequenceH const & seqH, TSequenceV const & seqV, int bandExtension)
{
    typedef typename Iterator<TSeedSet const, Standard>::Type TIterator;
    typedef typename Value<TSeedSet>::Type TSeed;

    SEQAN_ASSERT_GT_MSG(length(seedSet), 0u, "SeedSet is empty!");

    TIterator it = end(seedSet, Standard());
    --it;

    while (it != iterBegin)
    {
        TSeed seed = value(--it);
        if (_isCrossingEndHorizontal(endPositionH(seed), bandExtension, seqH))
            continue;
        else if (_isCrossingEndVertical(endPositionV(seed), bandExtension, seqV))
            continue;
        else  // Found seed which is not crossing the end.
            return it;
    }
    return it;
}

// ----------------------------------------------------------------------------
// Function _printTraceSegments()
// ----------------------------------------------------------------------------

template <typename T>
void _printTraceSegments(T const & globalTraceSet)
{
    for (unsigned i = 0; i < length(globalTraceSet); ++i)
    {
        std::cerr << "Tracback No " << i + 1 << std::endl;
        for (unsigned j = length(globalTraceSet[i]); j > 0; --j)
        {
            std::cerr << globalTraceSet[i][j-1] << "  ";
        }
        std::cerr << std::endl;
    }
}

// ----------------------------------------------------------------------------
// Function _initiaizeBeginningOfBandedChain()
// ----------------------------------------------------------------------------

// Handles the initialization of cells for the very first dp matrix.
template <typename TScoutState, typename TSizeH, typename TSizeV, typename TScoreScheme, typename TFreeEndGaps,
          typename TDPMatrixLocation, typename TGapCosts, typename TTraceback>
inline void
_initiaizeBeginningOfBandedChain(TScoutState & scoutState,
                                 TSizeH sizeH,
                                 TSizeV sizeV,
                                 TScoreScheme const & scoreScheme,
                                 DPProfile_<BandedChainAlignment_<TFreeEndGaps, TDPMatrixLocation>, TGapCosts, TTraceback> const &)
{
    typedef typename TScoutState::TInitCell TInitCell;
    typedef typename Value<TInitCell, 3>::Type TDPCell;
    typedef DPProfile_<BandedChainAlignment_<TFreeEndGaps, TDPMatrixLocation>, TGapCosts, TTraceback> TDPProfile;

    TDPCell dpInitCellHorizontal;
    _computeScore(dpInitCellHorizontal, TDPCell(), TDPCell(), TDPCell(), Nothing(), Nothing(), Nothing(),
                  RecursionDirectionZero(), TDPProfile());
    scoutState._nextInitializationCells.insert(TInitCell(0,0, dpInitCellHorizontal));

    for (TSizeH activeColumn = 1; activeColumn < sizeH; ++activeColumn)
    {
        TDPCell prevCell = dpInitCellHorizontal;
        if (IsFreeEndGap_<TFreeEndGaps, DPFirstRow>::VALUE)
            _computeScore(dpInitCellHorizontal, TDPCell(), TDPCell(), TDPCell(), Nothing(), Nothing(), scoreScheme,
                  RecursionDirectionZero(), TDPProfile());
        else
            _computeScore(dpInitCellHorizontal, TDPCell(), prevCell, TDPCell(), Nothing(), Nothing(), scoreScheme,
                  RecursionDirectionHorizontal(), TDPProfile());
        scoutState._nextInitializationCells.insert(TInitCell(activeColumn, 0, dpInitCellHorizontal));
    }

    TDPCell dpInitCellVertical;
    _computeScore(dpInitCellVertical, TDPCell(), TDPCell(), TDPCell(), Nothing(), Nothing(), Nothing(),
                  RecursionDirectionZero(), TDPProfile());
    for(TSizeV activeRow = 1; activeRow < sizeV; ++activeRow)
    {
        TDPCell prevCell = dpInitCellVertical;
        if (IsFreeEndGap_<TFreeEndGaps, DPFirstColumn>::VALUE)
            _computeScore(dpInitCellVertical, TDPCell(), TDPCell(), TDPCell(), Nothing(), Nothing(), scoreScheme,
                  RecursionDirectionZero(), TDPProfile());
        else
            _computeScore(dpInitCellVertical, TDPCell(), TDPCell(), prevCell, Nothing(), Nothing(), scoreScheme,
                  RecursionDirectionVertical(), TDPProfile());
        scoutState._nextInitializationCells.insert(TInitCell(0, activeRow, dpInitCellVertical));
    }

}

// ----------------------------------------------------------------------------
// Function _initializeBandedChain()
// ----------------------------------------------------------------------------

// Initializes the banded chain alignment. It computes the first band or the
// first gap-area followed by the first anchor.
template <typename TTraceSet, typename TScoutState, typename TSeed, typename TSeqH, typename TSeqV,
          typename TScoreScheme, typename TFreeEndGaps, typename TDPMatrixLocation, typename TGaps,
          typename TTracebackConfig>
inline typename Value<TScoreScheme>::Type
_initializeBandedChain(TTraceSet & globalTraceSet,
                       TScoutState & scoutState,
                       TSeed const & seed,
                       unsigned bandExtension,
                       TSeqH const & seqH,
                       TSeqV const & seqV,
                       TScoreScheme const & scoreSchemeAnchor,
                       TScoreScheme const & scoreSchemeGap,
                       DPProfile_<BandedChainAlignment_<TFreeEndGaps, TDPMatrixLocation>, TGaps, TracebackOn<TTracebackConfig> > const & dpProfile)
{
    typedef typename Infix<TSeqH const>::Type TInfixH;
    typedef typename Infix<TSeqV const>::Type TInfixV;
    typedef typename Position<TSeed const>::Type TSeedPosition;
    typedef typename MakeSigned<TSeedPosition>::Type TSignedPosition;
    typedef typename Size<TSeed const>::Type TSeedSize;
    typedef typename MakeSigned<TSeedSize>::Type TSignedSize;
    //typedef typename Value<TTraceSet>::Type TTraceArray;
    typedef typename Value<TScoreScheme>::Type TScore;
    typedef Pair<TSeedPosition, TSeedPosition> TGridPoint;

    // TODO(rmaerker): Change _horizontalBandShift to _horizontalBandShiftBeginPoint
    // TODO(rmaerker): Change _verticalBandShift to _verticalBandShiftBeginPoint
    TSeedSize horizontalBandShift = _horizontalBandShiftBeginPoint(seed);
    TSeedSize verticalBandShift = _verticalBandShiftBeginPoint(seed);

    TSeedSize horizontalNextGridOrigin = _max(0, static_cast<TSignedPosition>(beginPositionH(seed)) + 1 -
                                              static_cast<TSignedPosition>(bandExtension));
    TSeedSize verticalNextGridOrigin = _max(0, static_cast<TSignedPosition>(beginPositionV(seed)) + 1 -
                                            static_cast<TSignedPosition>(bandExtension));

    DPBandConfig<BandOn> band;
    // Determine coordinates of lower right corner of gap-area.
    setUpperDiagonal(band, static_cast<TSignedPosition>(_min(static_cast<TSignedSize>(length(seqH)),
        static_cast<TSignedPosition>(horizontalNextGridOrigin +  (bandExtension << 1) + horizontalBandShift +
        _max(0,static_cast<TSignedPosition>(bandExtension) - static_cast<TSignedPosition>(beginPositionV(seed)) -1)) +
        _min(0, static_cast<TSignedPosition>(beginPositionH(seed)) + 1 - static_cast<TSignedPosition>(bandExtension)))));

    setLowerDiagonal(band, -static_cast<TSignedPosition>(_min(static_cast<TSignedSize>(length(seqV)),
        static_cast<TSignedPosition>(verticalNextGridOrigin + (bandExtension << 1) + verticalBandShift +
        _max(0,static_cast<TSignedPosition>(bandExtension) - static_cast<TSignedPosition>(beginPositionH(seed)) -1)) +
        _min(0, static_cast<TSignedPosition>(beginPositionV(seed)) + 1 - static_cast<TSignedPosition>(bandExtension)))));

    TScore score = 0;
    TInfixH infixH;
    TInfixV infixV;

    // The first anchor does not cross the first row or column. Hence we have to fill the gap first.
    if (horizontalNextGridOrigin != 0u || verticalNextGridOrigin != 0u)
    {

        // INITIALIZATION
        _initiaizeBeginningOfBandedChain(scoutState, upperDiagonal(band) + 1, 1 - lowerDiagonal(band), scoreSchemeGap,
                                         dpProfile);
        // Define infixes covering the area which is to be computed.
        infixH = infix(seqH, 0, upperDiagonal(band));
        infixV = infix(seqV, 0, -lowerDiagonal(band));

        // Prepare the scout state for the next part of the algorithm.
        _reinitScoutState(scoutState, horizontalNextGridOrigin, verticalNextGridOrigin, 1 + upperDiagonal(band),
                          1 - lowerDiagonal(band), 1 + upperDiagonal(band) - horizontalNextGridOrigin,
                          1 - lowerDiagonal(band) - verticalNextGridOrigin);

        // Call the basic alignment function.
        score = _computeAlignment(globalTraceSet, scoutState, infixH, infixV, scoreSchemeGap, DPBandConfig<BandOff>(),
                                  DPProfile_<BandedChainAlignment_<TFreeEndGaps, BandedChainInitialDPMatrix>, TGaps,
                                             TracebackOn<TTracebackConfig> >());
    }
    else
    {
        _initiaizeBeginningOfBandedChain(scoutState, upperDiagonal(band), -lowerDiagonal(band), scoreSchemeAnchor,
                                         dpProfile);
    }

    TGridPoint gridBegin;
    // Determine coordinates of upper left corner of anchor.
    gridBegin.i1 = horizontalNextGridOrigin;  // This is either 0 if not used before or the origin of the current matrix.
    gridBegin.i2 = verticalNextGridOrigin;  // This is either 0 if not used before or the origin of the current matrix.

    // Check if the end point is at most as long as the sequences.
    TGridPoint gridEnd;
    gridEnd.i1 = _min(length(seqH), endPositionH(seed) + bandExtension);
    gridEnd.i2 = _min(length(seqV), endPositionV(seed) + bandExtension);

    // Define area covering the infixes.
    infixH = infix(seqH, gridBegin.i1, gridEnd.i1);
    infixV = infix(seqV, gridBegin.i2, gridEnd.i2);

    horizontalBandShift = _horizontalBandShiftEndPoint(seed);
    verticalBandShift = _verticalBandShiftEndPoint(seed);

    horizontalNextGridOrigin = _max(0, static_cast<TSignedPosition>(endPositionH(seed)) -
        static_cast<TSignedPosition>(bandExtension) - static_cast<TSignedPosition>(horizontalBandShift) -
        _max(0, static_cast<TSignedPosition>(endPositionV(seed) + bandExtension) -
        static_cast<TSignedSize>(length(seqV))) - static_cast<TSignedPosition>(gridBegin.i1));
    verticalNextGridOrigin = _max(0, static_cast<TSignedPosition>(endPositionV(seed)) -
        static_cast<TSignedPosition>(bandExtension) - static_cast<TSignedPosition>(verticalBandShift) -
        _max(0, static_cast<TSignedPosition>(endPositionH(seed) + bandExtension) -
        static_cast<TSignedSize>(length(seqH))) - static_cast<TSignedPosition>(gridBegin.i2));

    setUpperDiagonal(band, upperDiagonal(band) - static_cast<TSignedPosition>(gridBegin.i1));
    setLowerDiagonal(band, lowerDiagonal(band) + static_cast<TSignedPosition>(gridBegin.i2));

    // Calibrate the next vertical origin to the correct position within the band, only if it is a small band.
    if (static_cast<TSignedPosition>(length(infixV)) + lowerDiagonal(band) > upperDiagonal(band))
        verticalNextGridOrigin -= (length(infixV) + lowerDiagonal(band)) - upperDiagonal(band);

    _reinitScoutState(scoutState, horizontalNextGridOrigin, verticalNextGridOrigin, upperDiagonal(band) + 1,
                  1 - lowerDiagonal(band), length(infixH) - horizontalNextGridOrigin + 1,
                  length(infixV) - verticalNextGridOrigin + 1);

    // Prepare the scout state for the next part of the algorithm.
    TTraceSet localTraceSet;
    if (gridBegin.i1 == 0u && gridBegin.i2 == 0u)
    {
        if (gridEnd.i1 == length(seqH) && gridEnd.i2 == length(seqV))   // Standard global alignment problem.
        {
            resize(localTraceSet, 1);
            DPScoutState_<Default> noScout;
            score = _computeAlignment(localTraceSet[0], noScout, infixH, infixV, scoreSchemeGap, band,
                      DPProfile_<GlobalAlignment_<TFreeEndGaps>, TGaps, TracebackOn<TTracebackConfig> >());
        }
        else
            score = _computeAlignment(localTraceSet, scoutState, infixH, infixV, scoreSchemeGap, band,
                      DPProfile_<BandedChainAlignment_<TFreeEndGaps, BandedChainInitialDPMatrix>, TGaps, TracebackOn<TTracebackConfig> >());
    }
    else
    {
        if (gridEnd.i1 == length(seqH) && gridEnd.i2 == length(seqV))  // The anchor crosses the end of the matrix.
            score = _computeAlignment(localTraceSet, scoutState, infixH, infixV, scoreSchemeGap, band,
                              DPProfile_<BandedChainAlignment_<TFreeEndGaps, BandedChainFinalDPMatrix>, TGaps, TracebackOn<TTracebackConfig> >());
        else
            score = _computeAlignment(localTraceSet, scoutState, infixH, infixV, scoreSchemeGap, band, dpProfile);
    }

    if (gridBegin.i1 != 0u || gridBegin.i2 != 0u)
    {
        _adaptLocalTracesToGlobalGrid(localTraceSet, gridBegin);
        if (!empty(localTraceSet))
            _glueTracebacks(globalTraceSet, localTraceSet);
    }
    else
        globalTraceSet = localTraceSet;

    scoutState._horizontalNextGridOrigin += gridBegin.i1;
    scoutState._verticalNextGridOrigin += gridBegin.i2;

    if (static_cast<TSignedPosition>(length(infixV)) + lowerDiagonal(band) > upperDiagonal(band))
        scoutState._verticalNextGridOrigin += (static_cast<TSignedPosition>(length(infixV)) + lowerDiagonal(band)) -
                                               upperDiagonal(band);
    return score;
}

// ----------------------------------------------------------------------------
// Function _computeGapArea()
// ----------------------------------------------------------------------------

template <typename TTraceSet, typename TDPScoutState, typename TSeed, typename TSeqH, typename TSeqV,
          typename TScoreScheme, typename TFreeEndGaps, typename TDPMatrixLocation, typename TGaps, typename TTracebackSpec>
inline typename Value<TScoreScheme>::Type
_computeGapArea(TTraceSet & globalTraceSet,
                TDPScoutState & scoutState,
                TSeed const & currentSeed,
                unsigned bandExtension,
                TSeqH const & seqH,
                TSeqV const & seqV,
                TScoreScheme const & scoreScheme,
                DPProfile_<BandedChainAlignment_<TFreeEndGaps, TDPMatrixLocation>, TGaps, TracebackOn<TTracebackSpec> > const & dpProfile)
{
    typedef typename Infix<TSeqH const>::Type TInfixH;
    typedef typename Infix<TSeqV const>::Type TInfixV;
    typedef typename Position<TSeqH>::Type TPosH;
    typedef typename Position<TSeqV>::Type TPosV;
    typedef typename Position<TSeed const>::Type TPosition;
    typedef DPBandConfig<BandOff> TBand;
    typedef typename Value<TScoreScheme>::Type TScoreValue;
    typedef Pair<TPosition, TPosition> TGridPoint;

    TGridPoint gridBegin(scoutState._horizontalNextGridOrigin, scoutState._verticalNextGridOrigin);

    // TODO(rmaerker): Change _horizontalBandShift to _horizontalBandShiftBeginPoint
    // TODO(rmaerker): Change _verticalBandShift to _verticalBandShiftBeginPoint
    TGridPoint gridEnd(beginPositionH(currentSeed) + 1 + bandExtension + _horizontalBandShiftBeginPoint(currentSeed),
                       beginPositionV(currentSeed) + 1 + bandExtension + _verticalBandShiftBeginPoint(currentSeed));

    // Define infix area for alignment.
    TInfixH infixH = infix(seqH, gridBegin.i1, gridEnd.i1);
    TInfixV infixV = infix(seqV, gridBegin.i2, gridEnd.i2);

    // Relative begin positions of next grid.
    TPosH horizontalNextGridOrigin = beginPositionH(currentSeed) + 1 - bandExtension - gridBegin.i1;
    TPosV verticalNextGridOrigin = beginPositionV(currentSeed) + 1 - bandExtension - gridBegin.i2;

    // Prepare the scout state for the next part of the algorithm.
    _reinitScoutState(scoutState, horizontalNextGridOrigin, verticalNextGridOrigin, gridEnd.i1 - gridBegin.i1 + 1,
                      gridEnd.i2 - gridBegin.i2 + 1, gridEnd.i1 - gridBegin.i1 + 1 - horizontalNextGridOrigin,
                      gridEnd.i2 - gridBegin.i2 + 1 - verticalNextGridOrigin);

    // Compute the alignment.
    TTraceSet localTraceSet;
    TScoreValue score = _computeAlignment(localTraceSet, scoutState, infixH, infixV, scoreScheme, TBand(), dpProfile);

    // Adapt the local traces to match the positions of the  global grid.
    _adaptLocalTracesToGlobalGrid(localTraceSet, gridBegin);
    if (!empty(localTraceSet))
        _glueTracebacks(globalTraceSet, localTraceSet);

    scoutState._horizontalNextGridOrigin += gridBegin.i1;
    scoutState._verticalNextGridOrigin += gridBegin.i2;
    return score;
}

// ----------------------------------------------------------------------------
// Function _computeAnchorArea()
// ----------------------------------------------------------------------------

template <typename TTraceSet, typename TDPScoutState, typename TSeed, typename TSeqH, typename TSeqV,
          typename TScoreScheme, typename TFreeEndGaps, typename TDPMatrixLocation, typename TGaps, typename TTracebackSpec>
inline typename Value<TScoreScheme>::Type
_computeAnchorArea(TTraceSet & globalTraceSet,
                   TDPScoutState & scoutState,
                   TSeed const & currentSeed,
                   unsigned bandExtension,
                   TSeqH const & seqH,
                   TSeqV const & seqV,
                   TScoreScheme const & scoreScheme,
                   DPProfile_<BandedChainAlignment_<TFreeEndGaps, TDPMatrixLocation>, TGaps, TracebackOn<TTracebackSpec> > const & dpProfile)
{
    typedef typename Infix<TSeqH const>::Type TInfixH;
    typedef typename Infix<TSeqV const>::Type TInfixV;
    typedef typename Size<TSeed const>::Type TSize;
    typedef typename MakeSigned<TSize>::Type TSignedSize;
    typedef typename Position<TSeed const>::Type TPosition;
    typedef typename MakeSigned<TPosition>::Type TSignedPosition;
    //typedef DPBandConfig<BandOn> TBand;
    typedef typename Value<TScoreScheme>::Type TScore;
    typedef Pair<TPosition, TPosition> TGridPoint;

    TGridPoint gridBegin(scoutState._horizontalNextGridOrigin, scoutState._verticalNextGridOrigin);
    TGridPoint gridEnd(endPositionH(currentSeed) + bandExtension,
                       endPositionV(currentSeed) + bandExtension);

    // Define area covering the infixes of the current sub-matrix.
    TInfixH infixH = infix(seqH, gridBegin.i1, gridEnd.i1);
    TInfixV infixV = infix(seqV, gridBegin.i2, gridEnd.i2);

    TPosition horizontalBandShift = _horizontalBandShiftBeginPoint(currentSeed);
    TPosition verticalBandShift = _verticalBandShiftBeginPoint(currentSeed);

    // Computing the start point of the next grid that overlaps with the current anchor at the end.
    TPosition horizontalNextGridOrigin = endPositionH(currentSeed) - bandExtension - _horizontalBandShiftEndPoint(currentSeed) - gridBegin.i1;
    TPosition verticalNextGridOrigin = endPositionV(currentSeed) - bandExtension - _verticalBandShiftEndPoint(currentSeed) - gridBegin.i2;

    // Computing the start position of the band according to the start point.
    // TODO(rmaerker): Rename the function _verticalBandShift to _verticalBandShiftBeginPoint()
    // TODO(rmaerker): Rename the function _horiontalBandShift to _horizontalBandShiftBeginPoint()
    DPBandConfig<BandOn> band(-static_cast<TSignedPosition>((bandExtension << 1)) -static_cast<TSignedPosition>(verticalBandShift),
                         static_cast<TSignedPosition>((bandExtension << 1) + horizontalBandShift));

    TPosition relativeVerticalNextGridOrigin = verticalNextGridOrigin;
    // Calibrate the next vertical origin to the correct position within the band.
    if (static_cast<TSignedSize>(length(infixV)) + lowerDiagonal(band) > upperDiagonal(band))
        relativeVerticalNextGridOrigin -= (length(infixV) + lowerDiagonal(band)) - upperDiagonal(band);

    _reinitScoutState(scoutState, horizontalNextGridOrigin, relativeVerticalNextGridOrigin, upperDiagonal(band) + 1,
                      1 - lowerDiagonal(band), length(infixH) - horizontalNextGridOrigin + 1,
                      length(infixV) - verticalNextGridOrigin + 1);

    // Compute the alignment.
    TTraceSet localTraceSet;
    TScore score = _computeAlignment(localTraceSet, scoutState, infixH, infixV, scoreScheme, band, dpProfile);

    // Adapt the local traces to match the positions of the  global grid.
    _adaptLocalTracesToGlobalGrid(localTraceSet, gridBegin);
    if (!empty(localTraceSet))
        _glueTracebacks(globalTraceSet, localTraceSet);

    scoutState._horizontalNextGridOrigin += gridBegin.i1;
    scoutState._verticalNextGridOrigin += gridBegin.i2;
    if (static_cast<TSignedSize>(length(infixV)) + lowerDiagonal(band) > upperDiagonal(band))
        scoutState._verticalNextGridOrigin += (length(infixV) + lowerDiagonal(band)) - upperDiagonal(band);
    return score;
}

// ----------------------------------------------------------------------------
// Function _finishBandedChain()
// ----------------------------------------------------------------------------

template <typename TTraceSet, typename TDPScoutState, typename TSeed, typename TSeqH, typename TSeqV,
          typename TScoreScheme, typename TFreeEndGaps, typename TDPMatrixLocation, typename TGaps, typename TTracebackConfig>
inline typename Value<TScoreScheme>::Type
_finishBandedChain(TTraceSet & globalTraceSet,
                   TDPScoutState & dpScoutState,
                   TSeed const & currentSeed,
                   unsigned bandExtension,
                   TSeqH const & seqH,
                   TSeqV const & seqV,
                   TScoreScheme const & scoreSchemeAnchor,
                   TScoreScheme const & scoreSchemeGap,
                   DPProfile_<BandedChainAlignment_<TFreeEndGaps, TDPMatrixLocation>, TGaps, TracebackOn<TTracebackConfig> > const & dpProfile)
{
    //typedef typename Position<TSeqH const>::Type TPosH;
    //typedef typename Position<TSeqV const>::Type TPosV;
    typedef typename Position<TSeed const>::Type TPosition;
    typedef typename MakeSigned<TPosition>::Type TSignedPosition;
    typedef typename Size<TSeed const>::Type TSize;
    typedef typename MakeSigned<TSize>::Type TSignedSize;
    typedef typename Infix<TSeqH const>::Type TInfixH;
    typedef typename Infix<TSeqV const>::Type TInfixV;
    typedef typename Value<TScoreScheme>::Type TScoreValue;
    typedef Pair<TPosition, TPosition> TGridPoint;

    // Part A: Compute last connecting rectangle between last two anchors
    TGridPoint gridBegin(dpScoutState._horizontalNextGridOrigin, dpScoutState._verticalNextGridOrigin);

    // TODO(rmaerker): Change _horizontalBandShift to _horizontalBandShiftBeginPoint
    // TODO(rmaerker): Change _verticalBandShift to _verticalBandShiftBeginPoint
    TPosition horizontalBandShift = _horizontalBandShiftBeginPoint(currentSeed);
    TPosition verticalBandShift = _verticalBandShiftBeginPoint(currentSeed);

    TGridPoint gridEnd(_min(length(seqH), beginPositionH(currentSeed) + 1 + bandExtension + horizontalBandShift),
                       _min(length(seqV), beginPositionV(currentSeed) + 1 + bandExtension + verticalBandShift));

    // Define infix area for alignment.
    TInfixH infixH = infix(seqH, gridBegin.i1, gridEnd.i1);
    TInfixV infixV = infix(seqV, gridBegin.i2, gridEnd.i2);

    // Relative begin positions of next grid.
    TPosition horizontalNextGridOrigin = _max(0, static_cast<TSignedPosition>(beginPositionH(currentSeed)) + 1 -
        static_cast<TSignedPosition>(bandExtension) - static_cast<TSignedPosition>(gridBegin.i1));
    TPosition verticalNextGridOrigin = _max(0, static_cast<TSignedPosition>(beginPositionV(currentSeed)) + 1 -
        static_cast<TSignedPosition>(bandExtension) - static_cast<TSignedPosition>(gridBegin.i2));

    // Prepare the scout state for the next part of the algorithm.
    _reinitScoutState(dpScoutState, horizontalNextGridOrigin, verticalNextGridOrigin, length(infixH) + 1,
                      length(infixV) + 1, length(infixH) - horizontalNextGridOrigin + 1,
                      length(infixV) - verticalNextGridOrigin + 1);

    // Compute the alignment.
    TTraceSet localTraceSet;
    _computeAlignment(localTraceSet, dpScoutState, infixH, infixV, scoreSchemeGap, DPBandConfig<BandOff>(), dpProfile);

    // Adapt the local traces to match the positions of the  global grid.
    _adaptLocalTracesToGlobalGrid(localTraceSet, gridBegin);
    if (!empty(localTraceSet))
        _glueTracebacks(globalTraceSet, localTraceSet);

    // Part B: Compute the last anchor

    // Get the absolute positions of the beginning point of the band.
    gridBegin.i1 += horizontalNextGridOrigin;
    gridBegin.i2 += verticalNextGridOrigin;

    // Get the absolute positions of the end point of the band.
    gridEnd.i1 = _min(length(seqH), endPositionH(currentSeed) + bandExtension);
    gridEnd.i2 = _min(length(seqV), endPositionV(currentSeed) + bandExtension);

    infixH = infix(seqH, gridBegin.i1, gridEnd.i1);
    infixV = infix(seqV, gridBegin.i2, gridEnd.i2);

    // The last anchor is crossing the end of the global matrix in both directions.
    if(gridEnd.i1 == length(seqH) && gridEnd.i2 == length(seqV))
    {
        DPBandConfig<BandOn> band(-static_cast<TSignedPosition>(gridEnd.i2 - gridBegin.i2), static_cast<TSignedPosition>(gridEnd.i1 - gridBegin.i1));
        _reinitScoutState(dpScoutState, 0, 0, upperDiagonal(band)+ 1, 1 - lowerDiagonal(band),
                          upperDiagonal(band)+ 1, 1 - lowerDiagonal(band));
            // TODO(rmaerker): Should we not set the nextGridOrigin to 0 as it is the case for the last rectangle? We want to compute the last full path.
            // Compute the last anchor which crosses the end of the global grid.
        TScoreValue score = _computeAlignment(localTraceSet, dpScoutState, infixH, infixV, scoreSchemeAnchor, band,
                                  DPProfile_<BandedChainAlignment_<TFreeEndGaps, BandedChainFinalDPMatrix>, TGaps, TracebackOn<TTracebackConfig> >());

        _adaptLocalTracesToGlobalGrid(localTraceSet, gridBegin);
        if (!empty(localTraceSet))
            _glueTracebacks(globalTraceSet, localTraceSet);
        return score;
    }

    DPBandConfig<BandOn> band(-static_cast<TSignedPosition>((bandExtension << 1)) -static_cast<TSignedPosition>(verticalBandShift),
                         static_cast<TSignedPosition>((bandExtension << 1) + horizontalBandShift));

    // Compute the last anchor and the rectangle to fill the gap between the end of the global matrix and the last anchor.
    horizontalNextGridOrigin = _max(0,
        static_cast<TSignedPosition>(endPositionH(currentSeed)) - static_cast<TSignedPosition>(bandExtension) -
        static_cast<TSignedPosition>(_horizontalBandShiftEndPoint(currentSeed)) - static_cast<TSignedPosition>(gridBegin.i1) -
        _max(0, static_cast<TSignedPosition>(endPositionV(currentSeed) + bandExtension) - static_cast<TSignedSize>(length(seqV))));

    verticalNextGridOrigin = _max(0,
        static_cast<TSignedPosition>(endPositionV(currentSeed)) - static_cast<TSignedPosition>(bandExtension) -
        static_cast<TSignedPosition>(_verticalBandShiftEndPoint(currentSeed)) - static_cast<TSignedPosition>(gridBegin.i2) -
        _max(0, static_cast<TSignedPosition>(endPositionH(currentSeed) + bandExtension) - static_cast<TSignedSize>(length(seqH))));

    // Calibrate the next vertical origin to the correct position within the band.
    // TODO(rmaerker): Does this also apply to the case, when the anchor is crossing the end of the matrix?
    if (static_cast<TSignedSize>(length(infixV)) + lowerDiagonal(band) > upperDiagonal(band))
        verticalNextGridOrigin -= (length(infixV) + lowerDiagonal(band)) - upperDiagonal(band);

   _reinitScoutState(dpScoutState, horizontalNextGridOrigin, verticalNextGridOrigin, upperDiagonal(band) + 1,
                  1 - lowerDiagonal(band), length(infixH) - horizontalNextGridOrigin + 1,
                  length(infixV) - verticalNextGridOrigin + 1);

    clear(localTraceSet);
    // Compute the last anchor.
    TScoreValue score = _computeAlignment(localTraceSet, dpScoutState, infixH, infixV, scoreSchemeAnchor, band, dpProfile);
    _adaptLocalTracesToGlobalGrid(localTraceSet, gridBegin);
    if (!empty(localTraceSet))
        _glueTracebacks(globalTraceSet, localTraceSet);

    // Close the gap between last anchor and the end of the global grid.
    gridBegin.i1 += horizontalNextGridOrigin;
    if (static_cast<TSignedSize>(length(infixV)) + lowerDiagonal(band) > upperDiagonal(band))
        verticalNextGridOrigin += (length(infixV) + lowerDiagonal(band)) - upperDiagonal(band);
    gridBegin.i2 += verticalNextGridOrigin;

    _reinitScoutState(dpScoutState, 0, 0, length(seqH) - gridBegin.i1 + 1, length(seqV) - gridBegin.i2 + 1,
                      length(seqH) - gridBegin.i1 + 1, length(seqV) - gridBegin.i2 + 1);

    // Compute the alignment.
    clear(localTraceSet);
    score = _computeAlignment(localTraceSet, dpScoutState, suffix(seqH, gridBegin.i1), suffix(seqV,gridBegin.i2),
                              scoreSchemeGap, DPBandConfig<BandOff>(), DPProfile_<BandedChainAlignment_<TFreeEndGaps,
                              BandedChainFinalDPMatrix>, TGaps, TracebackOn<TTracebackConfig> >());

    _adaptLocalTracesToGlobalGrid(localTraceSet, gridBegin);
    if (!empty(localTraceSet))
        _glueTracebacks(globalTraceSet, localTraceSet);
    return score;
}

// ----------------------------------------------------------------------------
// Function _computeAlignment()                                    [DP-Wrapper]
// ----------------------------------------------------------------------------

template <typename TGapScheme, typename TTraceTarget, typename TScoutState, typename TSequenceH, typename TSequenceV,
          typename TScoreScheme, typename TBandSwitch, typename TAlignmentAlgorithm, typename TTraceFlag>
inline typename Value<TScoreScheme>::Type
_computeAlignment(TTraceTarget & traceSegments,
                  TScoutState & scoutState,
                  TSequenceH const & seqH,
                  TSequenceV const & seqV,
                  TScoreScheme const & scoreScheme,
                  DPBandConfig<TBandSwitch> const & band,
                  DPProfile_<TAlignmentAlgorithm, TGapScheme, TTraceFlag> const & dpProfile)
{
    typedef typename Value<TScoreScheme>::Type TScoreValue;
    typedef DPContext<TScoreValue, TGapScheme> TDPContext;
    TDPContext dpContext;
    return _computeAlignment(dpContext, traceSegments, scoutState, seqH, seqV, scoreScheme, band, dpProfile);
}

// ----------------------------------------------------------------------------
// Function _computeAlignment()                          [BandedChainAlignment]
// ----------------------------------------------------------------------------

// Given a monotonic non-decreasing seed chain the following algorithms computes a
// global alignment around this rough global map, while connecting two anchors  by filling
// the gap between them using a standard dp algorithm.
template <typename TTraceSet, typename TSeedSet, typename TSequenceH, typename TSequenceV, typename TScoreValue,
          typename TScoreSpecAnchor, typename TScoreSpecGap, typename TFreeEndGaps, typename TDPMatrixLocation,
          typename TGapSpec, typename TTracebackConfig>
inline TScoreValue
_computeAlignment(TTraceSet & globalTraceSet,
                  TSeedSet const & seedSet,
                  TSequenceH const & seqH,
                  TSequenceV const & seqV,
                  Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchor,
                  Score<TScoreValue, TScoreSpecGap> const & scoreSchemeGap,
                  unsigned bandExtension,
                  DPProfile_<BandedChainAlignment_<TFreeEndGaps, TDPMatrixLocation>, TGapSpec,
                             TracebackOn<TTracebackConfig> > const & profile)
{
    //typedef DPProfile_<BandedChainAlignment_<TFreeEndGaps, TDPMatrixLocation>, TGapSpec, TracebackOn<TTracebackConfig> > TAlignmentProfile;

    typedef typename Position<TSequenceH>::Type TPosH;
    typedef typename Position<TSequenceV>::Type TPosV;
    typedef typename Iterator<TSeedSet const, Standard>::Type TSeedSetIterator;

    typedef DPCell_<TScoreValue, TGapSpec> TDPCell;
    typedef DPScoutState_<BandedChainAlignmentScoutState<TDPCell> > TScoutState;

    // handle cases of empty sequences
    SEQAN_ASSERT_MSG(!empty(seqH), "Cannot compute alignment on empty sequence (horizontal sequence)");
    SEQAN_ASSERT_MSG(!empty(seqV), "Cannot compute alignment on empty sequence (vertical sequence)");

    // TODO(rmaerker): At the moment the band must have at least a minimalBandwidth of size three. Maybe we change that.
    SEQAN_ASSERT_GEQ(bandExtension, 1u);

    // Handle case of empty seed set.
    if (length(seedSet) < 1)
        return MinValue<TScoreValue>::VALUE;

    // Find the first anchor that is not covered by the region between the beginning of the matrix and the next anchor
    // considering the minbandwidh parameter.
    TSeedSetIterator it = _findFirstAnchor(seedSet, bandExtension);
    // Find the last anchor that is not covered by the region between the previous anchor and the end of the matrix.
    TSeedSetIterator itEnd = _findLastAnchor(it, seedSet, seqH, seqV, bandExtension);

    SEQAN_ASSERT(itEnd != static_cast<TSeedSetIterator>(end(seedSet, Standard())));
    // The scout state stores the current state of the dp scout and is used to store the
    // initialization values for the next intersecting grid.
    TScoutState scoutState;

    // INITIALIZATION
    // Initialize the banded chain alignment to compute the first gap.
    TScoreValue score = _initializeBandedChain(globalTraceSet, scoutState, value(it), bandExtension, seqH, seqV,
                                               scoreSchemeAnchor, scoreSchemeGap, profile);

    // We need to close the band already if there is no other seed to consider.
    if (length(seedSet) == 1 || (it == itEnd && _isLastSeed(itEnd, seedSet)))
    {
        // There is nothing to close.
        if (endPositionH(value(it)) + bandExtension < length(seqH) ||
            endPositionV(value(it)) + bandExtension < length(seqV))
        {
            TPosH gridBeginH = scoutState._horizontalNextGridOrigin;
            TPosV gridBeginV = scoutState._verticalNextGridOrigin;
            _reinitScoutState(scoutState, 0, 0, length(seqH) + 1 - gridBeginH, length(seqV) + 1 - gridBeginV,
                              length(seqH) + 1 - gridBeginH, length(seqV) + 1 - gridBeginV);

            // compute the alignment
            TTraceSet localTraceSet;
            score = _computeAlignment(localTraceSet, scoutState, suffix(seqH, gridBeginH), suffix(seqV, gridBeginV),
                              scoreSchemeGap, DPBandConfig<BandOff>(), DPProfile_<BandedChainAlignment_<TFreeEndGaps,
                              BandedChainFinalDPMatrix>, TGapSpec, TracebackOn<TTracebackConfig> >());
            _adaptLocalTracesToGlobalGrid(localTraceSet, Pair<TPosH, TPosV>(gridBeginH, gridBeginV));
            _glueTracebacks(globalTraceSet, localTraceSet);
        }
        return score;
    }

    // MAIN: Process all inner seeds.
    while (it != itEnd)
    {
        SEQAN_ASSERT(_checkColinearity(it));
        // Process the next gap area that connects to anchors.
        _computeGapArea(globalTraceSet, scoutState, value(++it), bandExtension, seqH, seqV, scoreSchemeGap,
                        profile);
        // Process the folowing anchor.
        _computeAnchorArea(globalTraceSet, scoutState, value(it), bandExtension, seqH, seqV, scoreSchemeAnchor, profile);
    }
    SEQAN_ASSERT(_checkColinearity(it));
    // Finish the banded chain alignment while computing the closing gap.
    score = _finishBandedChain(globalTraceSet, scoutState, value(++it), bandExtension, seqH, seqV, scoreSchemeAnchor,
                               scoreSchemeGap, profile);
    return score;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_SEEDS_BANDED_CHAIN_ALIGNMENT_IMPL_H_
