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
// The dp scout is a structure that stores the current maximal score and its
// host position in the underlying dp-matrix.
// This class can be overloaded to implement different behaviors of tracking
// the maximal score, e.g., for the split breakpoint computation.
// ==========================================================================
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// The terminator specialization of dp scout offers the possibility to have
// dp generation stop, if specified criteria are met.
// To do this, define HasTerminationCriterium_<> for your algorithm and
// implement a DPScoutState for your terminator specialization. In your
// overloaded _scoutBestScore() or _computeCell() you can call
// terminateScout() on your Scout to have DP-generation stop.
// see dp_scout_xdrop.h for an example.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_TEST_ALIGNMENT_DP_SCOUT_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_TEST_ALIGNMENT_DP_SCOUT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Terminator_
// ----------------------------------------------------------------------------

template <typename TSpec = void>
struct Terminator_;

// ----------------------------------------------------------------------------
// Class DPScoutState_
// ----------------------------------------------------------------------------

template <typename TSpec>
class DPScoutState_;

template <>
class DPScoutState_<Default> : public Nothing  // empty member optimization
{};

// ----------------------------------------------------------------------------
// Class DPScout_
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TSpec>
class DPScout_;

// The default implementation of the dp scout simply stores one maximum
// and its corresponding position.
//
// The state must be a Nothing and is left untouched and unused.
template <typename TDPCell, typename TSpec>
class DPScout_
{
public:
    typedef typename Value<TDPCell>::Type TScoreValue;
//    TScoreValue _maxScore;          // The maximal score.
    TDPCell _maxScore;
    unsigned int _maxHostPosition;  // The corresponding host position within the underlying dp-matrix.

    DPScout_() : _maxScore(), _maxHostPosition(0) {}

    DPScout_(DPScoutState_<Default> const & /*state*/) :
        _maxScore(), _maxHostPosition(0) {}

    DPScout_(DPScout_ const & other) :
        _maxScore(other._maxScore), _maxHostPosition(other._maxHostPosition) {}

    DPScout_ & operator=(DPScout_ const & other)
    {
        if (this != &other)
        {
            _maxScore = other._maxScore;
            _maxHostPosition = other._maxHostPosition;
        }
        return *this;
    }

};

// Terminator_ Specialization
template <typename TDPCell, typename TSpec>
class DPScout_<TDPCell, Terminator_<TSpec> >
    : public DPScout_<TDPCell, Default>
{
public:
    typedef DPScout_<TDPCell, Default>  TParent;
    bool terminationCriteriumMet;
    DPScoutState_<Terminator_<TSpec> > * state;

    DPScout_()
        : TParent(),
          terminationCriteriumMet(false),
          state(0)
    {}

    DPScout_(DPScoutState_<Terminator_<TSpec> > & state)
        : TParent(),
          terminationCriteriumMet(false),
          state(&state)
    {}

    DPScout_(DPScout_ const & other)
        : TParent(static_cast<TParent const &>(other)),
          terminationCriteriumMet(other.terminationCriteriumMet),
          state(other.state)
    {}

    DPScout_ & operator=(DPScout_ const & other)
    {
        if (this != &other)
        {
            *static_cast<TParent*>(this) = other;
            terminationCriteriumMet = other.terminationCriteriumMet;
            state = other.state;
        }
        return *this;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ScoutSpecForAlignmentAlgorithm_
// ----------------------------------------------------------------------------

// Given an alignment algorithm tag such as GlobalAlignment_ or LocalAlignment_, returns the specialization tag for the
// corresponding DPScout_ specialization.

template <typename TAlignmentAlgorithm>
struct ScoutSpecForAlignmentAlgorithm_
{
    typedef If<HasTerminationCriterium_<TAlignmentAlgorithm>,
               Terminator_<>,
               Default> Type;
};

// ----------------------------------------------------------------------------
// Metafunction ScoutStateSpecForScout_
// ----------------------------------------------------------------------------

// Given an dp scout this meta-function returns the appropriate specialization for the scout state.

template <typename TScout>
struct ScoutStateSpecForScout_
{
    typedef Default Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _scoutBestScore()
// ----------------------------------------------------------------------------

// Tracks the new score, if it is the new maximum.
template <typename TDPCell, typename TSpec, typename TTraceMatrixNavigator,
          typename TIsLastColumn, typename TIsLastRow>
inline void
_scoutBestScore(DPScout_<TDPCell, TSpec> & dpScout,
                TDPCell const & activeCell,
                TTraceMatrixNavigator const & navigator,
                TIsLastColumn const & /**/,
                TIsLastRow const & /**/)
{

    if (_scoreOfCell(activeCell) > _scoreOfCell(dpScout._maxScore))
    {
        dpScout._maxScore = activeCell;
        dpScout._maxHostPosition = position(navigator);
    }
}

// TODO(rmaerker): Why is this needed?
template <typename TDPCell, typename TSpec, typename TTraceMatrixNavigator, typename TIsLastColumn>
inline void
_scoutBestScore(DPScout_<TDPCell, TSpec> & dpScout,
                TDPCell const & activeCell,
                TTraceMatrixNavigator const & navigator,
                TIsLastColumn const & /**/)
{
    return _scoutBestScore(dpScout, activeCell, navigator, TIsLastColumn(),
                           False());
}

// TODO(rmaerker): Why is this needed?
template <typename TDPCell, typename TSpec, typename TTraceMatrixNavigator>
inline void
_scoutBestScore(DPScout_<TDPCell, TSpec> & dpScout,
                TDPCell const & activeCell,
                TTraceMatrixNavigator const & navigator)
{
    return _scoutBestScore(dpScout, activeCell, navigator, False(), False());
}

// ----------------------------------------------------------------------------
// Function maxScore()
// ----------------------------------------------------------------------------

// Returns the current maximal score.
template <typename TDPCell, typename TScoutSpec>
inline typename Value<TDPCell>::Type const
maxScore(DPScout_<TDPCell, TScoutSpec> const & dpScout)
{
    return _scoreOfCell(dpScout._maxScore);
}

// ----------------------------------------------------------------------------
// Function maxHostPosition()
// ----------------------------------------------------------------------------

// Returns the host position that holds the current maximum score.
template <typename TDPCell, typename TScoutSpec>
inline unsigned int
maxHostPosition(DPScout_<TDPCell, TScoutSpec> const & dpScout)
{
    return dpScout._maxHostPosition;
}

// ----------------------------------------------------------------------------
// Function _terminationCriteriumIsMet()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TSpec>
inline bool
_terminationCriteriumIsMet(DPScout_<TDPCell, Terminator_<TSpec> > const & scout)
{
    return scout.terminationCriteriumMet;
}

// ----------------------------------------------------------------------------
// Function terminateScout()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TSpec>
inline void
terminateScout(DPScout_<TDPCell, Terminator_<TSpec> > & scout)
{
    scout.terminationCriteriumMet = true;
}



}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_TEST_ALIGNMENT_DP_SCOUT_H_
