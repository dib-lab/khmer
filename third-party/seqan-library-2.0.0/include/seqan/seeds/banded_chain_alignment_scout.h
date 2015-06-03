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
// Implements the DPScout for the banded chain alignment algorithm.
// The ScoutState is used to keep track of the initialization values for
// the next matrix.
// ==========================================================================

#ifndef INCLUDE_SEQAN_SEEDS_BANDED_CHAIN_ALIGNMENT_SCOUT_H_
#define INCLUDE_SEQAN_SEEDS_BANDED_CHAIN_ALIGNMENT_SCOUT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag BandedChainAlignmentScout
// ----------------------------------------------------------------------------

struct BandedChainAlignmentScout_;
typedef Tag<BandedChainAlignmentScout_> BandedChainAlignmentScout;

// ----------------------------------------------------------------------------
// Class BandedChainAlignmentScoutState
// ----------------------------------------------------------------------------

// Used to determine the scout state of the BandedChainAlignmentScout.
template <typename TSpec>
struct BandedChainAlignmentScoutState;

// ----------------------------------------------------------------------------
// Class DPScoutState_
// ----------------------------------------------------------------------------

// Stores the state of the algorithm.
// It keeps track of the initialization values for horziontal and vertical
// direction of the current matrix and the next matrix that intersects with
// the current active matrix.
template <typename TDPCell>
class DPScoutState_<BandedChainAlignmentScoutState<TDPCell> >
{
public:
    typedef Triple<unsigned, unsigned, TDPCell> TInitCell;
    typedef std::set<TInitCell> TInitializationCellSet;

    unsigned int _horizontalNextGridOrigin;
    unsigned int _verticalNextGridOrigin;

    String<TDPCell> _horizontalInitCurrentMatrix;
    String<TDPCell> _verticalInitCurrentMatrix;
    String<TDPCell> _horizontalInitNextMatrix;
    String<TDPCell> _verticalInitNextMatrix;
    TInitializationCellSet _nextInitializationCells;

    DPScoutState_() : _horizontalNextGridOrigin(0u), _verticalNextGridOrigin(0u)
    {}

};

// ----------------------------------------------------------------------------
// Class DPScout_
// ----------------------------------------------------------------------------

template<typename TDPCell>
class DPScout_<TDPCell, BandedChainAlignmentScout>
{
public:
    typedef typename ScoutStateSpecForScout_<DPScout_>::Type TDPScoutStateSpec;
    typedef DPScoutState_<TDPScoutStateSpec> TScoutState;
    typedef String<unsigned int> TMaxHostPositionString;
    typedef typename Value<TDPCell>::Type TScoreValue;

//    TScoreValue _maxScore;      // the maximal score detected
    TDPCell _maxScore;
    TScoutState * _dpScoutStatePtr;
    TMaxHostPositionString _maxHostPositions;    // the array containing all positions that have a maximal score

    DPScout_() : _maxScore(),
                _dpScoutStatePtr(0),
                _maxHostPositions()
    {}

    DPScout_(TScoutState & scoutState) :
                _maxScore(),
                _dpScoutStatePtr(&scoutState),
                _maxHostPositions()
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ScoutSpecForAlignmentAlgorithm_
// ----------------------------------------------------------------------------

template <typename TSpec, typename TDPMatrixLocation>
struct ScoutSpecForAlignmentAlgorithm_<BandedChainAlignment_<TSpec, TDPMatrixLocation> >
{
    typedef BandedChainAlignmentScout Type;
};

// ----------------------------------------------------------------------------
// Metafunction ScoutStateSpecForScout_
// ----------------------------------------------------------------------------

template <typename TDPCell>
struct ScoutStateSpecForScout_<DPScout_<TDPCell, BandedChainAlignmentScout> >
{
    typedef BandedChainAlignmentScoutState<TDPCell> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _rinitScout()
// ----------------------------------------------------------------------------

template <typename TDPCell>
inline void
_reinitScout(DPScout_<TDPCell, BandedChainAlignmentScout> & dpScout)
{
    //typedef typename Value<TDPCell>::Type TScoreValue;

    dpScout._maxScore = TDPCell();
    resize(dpScout._maxHostPositions, 1, 0);
}

// ----------------------------------------------------------------------------
// Function _reinitScoutState()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TPosH, typename TPosV, typename TSizeCurrInit, typename TSizeNextInit>
inline void
_reinitScoutState(DPScoutState_<BandedChainAlignmentScoutState<TDPCell> > & scoutState,
                  TPosH const & originNextMatrixH,
                  TPosV const & originNextMatrixV,
                  TSizeCurrInit const & sizeCurrMatrixInitH,
                  TSizeCurrInit const & sizeCurrMatrixInitV,
                  TSizeNextInit const & sizeNextMatrixInitH,
                  TSizeNextInit const & sizeNextMatrixInitV)
{
    typedef DPScoutState_<BandedChainAlignmentScoutState<TDPCell> > TDPScoutState;
    typedef typename TDPScoutState::TInitializationCellSet TInitCellSet;
    typedef typename TInitCellSet::iterator TInitCellSetIterator;

    scoutState._horizontalNextGridOrigin = originNextMatrixH;
    scoutState._verticalNextGridOrigin = originNextMatrixV;

    // initialize the current initialization values of the tracker
    arrayFill(begin(scoutState._horizontalInitCurrentMatrix), end(scoutState._horizontalInitCurrentMatrix), TDPCell());
    arrayFill(begin(scoutState._verticalInitCurrentMatrix), end(scoutState._verticalInitCurrentMatrix), TDPCell());
    arrayFill(begin(scoutState._horizontalInitNextMatrix), end(scoutState._horizontalInitNextMatrix), TDPCell());
    arrayFill(begin(scoutState._verticalInitNextMatrix), end(scoutState._verticalInitNextMatrix), TDPCell());

    // check if the value needs to be resized ... (can be longer but not smaller)
    if ((TSizeCurrInit) length(scoutState._horizontalInitCurrentMatrix) < sizeCurrMatrixInitH)
        resize(scoutState._horizontalInitCurrentMatrix, sizeCurrMatrixInitH, TDPCell());
    if ((TSizeCurrInit) length(scoutState._verticalInitCurrentMatrix) < sizeCurrMatrixInitV)
        resize(scoutState._verticalInitCurrentMatrix, sizeCurrMatrixInitV, TDPCell());

    if ((TSizeNextInit) length(scoutState._horizontalInitNextMatrix) < sizeNextMatrixInitH)
        resize(scoutState._horizontalInitNextMatrix, sizeNextMatrixInitH, TDPCell());
    if ((TSizeNextInit) length(scoutState._verticalInitNextMatrix) < sizeNextMatrixInitV)
        resize(scoutState._verticalInitNextMatrix, sizeNextMatrixInitV, TDPCell());

    // Parsing the set to get the values.
    for (TInitCellSetIterator it = scoutState._nextInitializationCells.begin();  // We need to plant the initialization values here. At the moment we only put the old values here.
        it != scoutState._nextInitializationCells.end(); ++it)
    {
//        std::cerr << "TInitCell == " << it->i1 << " " << it->i2 << " " << _scoreOfCell(it->i3) << std::endl;
        if (it->i1 == 0)
            scoutState._verticalInitCurrentMatrix[it->i2] = it->i3;
        if (it->i2 == 0)
            scoutState._horizontalInitCurrentMatrix[it->i1] = it->i3;

    }
}

// ----------------------------------------------------------------------------
// Function _setScoutState()
// ----------------------------------------------------------------------------

// Sets the state of the scout.
template <typename TDPCell, typename TStateSpec>
inline void
_setScoutState(DPScout_<TDPCell, BandedChainAlignmentScout> & dpScout,
               DPScoutState_<BandedChainAlignmentScoutState<TStateSpec> > & state)
{
    dpScout._dpScoutStatePtr = & state;
}

// ----------------------------------------------------------------------------
// Function _scoutBestScore()
// ----------------------------------------------------------------------------


template <typename TDPCell, typename TTraceMatrixNavigator>
inline void
_scoutBestScore(DPScout_<TDPCell, BandedChainAlignmentScout> & dpScout,
                TDPCell const & dpCell,
                TTraceMatrixNavigator const & navigator,
                bool isLastColumn,
                bool isLastRow,
                bool trackNextInitColumn,
                bool trackNextInitRow)
{
    // Store value for vertical initialization of next grid.
    if(trackNextInitColumn)
        dpScout._dpScoutStatePtr->_verticalInitNextMatrix[coordinate(navigator, +DPMatrixDimension_::VERTICAL) -
                                                       dpScout._dpScoutStatePtr->_verticalNextGridOrigin] = dpCell;
    // Store value for horizontal initialization of next grid.
    if (trackNextInitRow)
        dpScout._dpScoutStatePtr->_horizontalInitNextMatrix[coordinate(navigator, +DPMatrixDimension_::HORIZONTAL) -
                                                         dpScout._dpScoutStatePtr->_horizontalNextGridOrigin] = dpCell;

    // Now we have to track the optimal score from the last column or row.
    if (isLastColumn || isLastRow)
    {
        if (_scoreOfCell(dpCell) >= _scoreOfCell(dpScout._maxScore))
        {
            if (_scoreOfCell(dpCell) == _scoreOfCell(dpScout._maxScore))
                appendValue(dpScout._maxHostPositions, position(navigator));
            else
            {
                resize(dpScout._maxHostPositions, 1);
                dpScout._maxHostPositions[0] = position(navigator);
                dpScout._maxScore = dpCell;
            }
        }
    }
}

template <typename TDPCell, typename TTraceMatrixNavigator>
inline void
_scoutBestScore(DPScout_<TDPCell, BandedChainAlignmentScout> &,
                TDPCell const &,
                TTraceMatrixNavigator const &,
                bool isLastColumn = false,
                bool isLastRow = false)
{
    (void) isLastColumn;
    (void) isLastRow;
    //no-op
}

// Delegate.
template <typename TDPCell, typename TTraceMatrixNavigator, typename TIsLastColumn, typename TIsLastRow>
inline void
_scoutBestScore(DPScout_<TDPCell, BandedChainAlignmentScout> & dpScout,
                TDPCell const & activeCell,
                TTraceMatrixNavigator const & navigator,
                TIsLastColumn const & /**/,
                TIsLastRow const & /**/)
{
    _scoutBestScore(dpScout, activeCell, navigator, TIsLastColumn::VALUE, TIsLastRow::VALUE);
}

// ----------------------------------------------------------------------------
// Function maxScore()
// ----------------------------------------------------------------------------

template <typename TDPCell>
inline typename Value<TDPCell>::Type &
maxScore(DPScout_<TDPCell, BandedChainAlignmentScout> & scout)
{
    return _scoreOfCell(scout._maxScore);
}

template <typename TDPCell>
inline typename Value<TDPCell>::Type const &
maxScore(DPScout_<TDPCell, BandedChainAlignmentScout> const & scout)
{
    return _scoreOfCell(scout._maxScore);
}

// ----------------------------------------------------------------------------
// Function maxHostPositions()
// ----------------------------------------------------------------------------

template <typename TDPCell>
inline typename DPScout_<TDPCell, BandedChainAlignmentScout>::TMaxHostPositionString const &
maxHostPositions(DPScout_<TDPCell, BandedChainAlignmentScout> const & scout)
{
    return scout._maxHostPositions;
}

template <typename TDPCell>
inline typename DPScout_<TDPCell, BandedChainAlignmentScout>::TMaxHostPositionString &
maxHostPositions(DPScout_<TDPCell, BandedChainAlignmentScout> & scout)
{
    return scout._maxHostPositions;
}

// ----------------------------------------------------------------------------
// Function maxHostPosition()
// ----------------------------------------------------------------------------

template <typename TDPCell>
inline unsigned int
maxHostPosition(DPScout_<TDPCell, BandedChainAlignmentScout> const & scout)
{
    return scout._maxHostPositions[0];
}

// ----------------------------------------------------------------------------
// Function nextMatrixBeginH()
// ----------------------------------------------------------------------------

template <typename TDPCell>
inline unsigned int
_nextMatrixBeginH(DPScout_<TDPCell, BandedChainAlignmentScout> const & scout)
{
    return scout._posH;
}

// ----------------------------------------------------------------------------
// Function nextMatrixBeginV()
// ----------------------------------------------------------------------------

template <typename TDPCell>
inline unsigned int
_nextMatrixBeginV(DPScout_<TDPCell, BandedChainAlignmentScout> const & scout)
{
    return scout._posV;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_SEEDS_BANDED_CHAIN_ALIGNMENT_SCOUT_H_
