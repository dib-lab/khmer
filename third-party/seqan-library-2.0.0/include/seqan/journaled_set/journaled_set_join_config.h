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

#ifndef INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_JOIN_CONFIG_H_
#define INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_JOIN_CONFIG_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class JoinConfig
// ----------------------------------------------------------------------------

template <>
struct JoinConfig<GlobalAlign<JournaledManhatten> > {};

template <>
struct JoinConfig<GlobalAlign<JournaledCompact> >
{
    typedef Score<int, BiAffine> TScoringScheme;
    TScoringScheme _score;
    int  _lowerDiagonal;
    int  _upperDiagonal;
    bool _isBandSet;
    AlignConfig<true, false, false, true> _alignConfig;

    JoinConfig() : _score(0, -100000, 0, -12, -1, -25),
                   _lowerDiagonal(0),
                   _upperDiagonal(0),
                   _isBandSet(false),
                   _alignConfig()
    {}
};


// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function isBandSet()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool isBandSet(JoinConfig<TSpec> const & joinConfig)
{
    return joinConfig._isBandSet;
}

// ----------------------------------------------------------------------------
// Function lowerDiagonal()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline int lowerDiagonal(JoinConfig<TSpec> const & joinConfig)
{
    return joinConfig._lowerDiagonal;
}

// ----------------------------------------------------------------------------
// Function setLowerDiagonal()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void setLowerDiagonal(JoinConfig<TSpec> & joinConfig, int lowerDiagonal)
{
    joinConfig._isBandSet = true;
    joinConfig._lowerDiagonal = lowerDiagonal;
}

// ----------------------------------------------------------------------------
// Function setLowerDiagonal()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline int upperDiagonal(JoinConfig<TSpec> const & joinConfig)
{
    return joinConfig._upperDiagonal;
}

// ----------------------------------------------------------------------------
// Function setUpperDiagonal()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void setUpperDiagonal(JoinConfig<TSpec> & joinConfig, int upperDiagonal)
{
    joinConfig._isBandSet = true;
    joinConfig._upperDiagonal = upperDiagonal;
}


inline Score<int, BiAffine> &
scoringScheme(JoinConfig<GlobalAlign<JournaledCompact> > & joinConfig)
{
    return joinConfig._score;
}

inline Score<int, BiAffine> const &
scoringScheme(JoinConfig<GlobalAlign<JournaledCompact> > const & joinConfig)
{
    return joinConfig._score;
}

// ----------------------------------------------------------------------------
// Function setScoringScheme()
// ----------------------------------------------------------------------------

inline void
setScoringScheme(JoinConfig<GlobalAlign<JournaledCompact> > & joinConfig, Score<int, BiAffine> const & scoringScheme)
{
    joinConfig._score = scoringScheme;
    setScoreMismatch(joinConfig._score, -100000);   // Explicitly forbid mis matches in context with journaling.
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_JOIN_CONFIG_H_
