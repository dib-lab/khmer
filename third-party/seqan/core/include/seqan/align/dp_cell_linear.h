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
// Implements the DPCell for linear gap cost functions. This requires only
// one matrix, such that we only need to store one score value per matrix
// entry.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_CELL_LINEAR_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_CELL_LINEAR_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================


// ----------------------------------------------------------------------------
// Class DPCell                                                    [LinearGaps]
// ----------------------------------------------------------------------------

// The specialization for linear gap cost function.
// It solely stores the maximal score.
template <typename TScoreValue>
class DPCell_<TScoreValue, LinearGaps>
{
public:

    TScoreValue _score;

    // The default c'tor.
    DPCell_() :
        _score(DPCellDefaultInfinity<DPCell_>::VALUE)
    {}

    // The copy c'tor.
    DPCell_(DPCell_<TScoreValue, LinearGaps> const & other) :
        _score(other._score)
    {}

    // The assignment operator.
    DPCell_ &
    operator=(DPCell_<TScoreValue, LinearGaps> const & other)
    {
        if (this != &other)
            _score = other._score;
        return *this;
    }

    DPCell_ &
    operator=(TScoreValue const & score)
    {
        _score = score;
        return *this;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function operator<()
// ----------------------------------------------------------------------------

// Needed for banded chain alignment for the std::set.
template <typename TScoreValueLeft, typename TScoreValueRight>
inline bool operator<(DPCell_<TScoreValueLeft, LinearGaps> const & left,
                      DPCell_<TScoreValueRight, LinearGaps> const & right)
{
    return left._score < right._score;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_CELL_LINEAR_H_
