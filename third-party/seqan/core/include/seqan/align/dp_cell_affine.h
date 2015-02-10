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
// Implements the DPCell for affine gap functions. It stores the score
// values for the three matrices: diagonal, vertical and horizontal.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_CELL_AFFINE_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_CELL_AFFINE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class DPCell                                                    [AffineGaps]
// ----------------------------------------------------------------------------

// Specialization for the affine gap cost function.
// This class stores three values belonging to the matrix storing the current
// maximum, the matrix for horizontal gaps and the matrix for vertical gaps.
template <typename TScoreValue>
class DPCell_<TScoreValue, AffineGaps>
{
public:
    TScoreValue _score;
    TScoreValue _horizontalScore;
    TScoreValue _verticalScore;

    // The default c'tor.
    DPCell_() :
        _score(DPCellDefaultInfinity<DPCell_>::VALUE),
        _horizontalScore(DPCellDefaultInfinity<DPCell_>::VALUE),
        _verticalScore(DPCellDefaultInfinity<DPCell_>::VALUE)
    {}

    // The copy c'tor.
    DPCell_(DPCell_<TScoreValue, AffineGaps> const & other) :
        _score(other._score),
        _horizontalScore(other._horizontalScore),
        _verticalScore(other._verticalScore)
    {}

    // The assignment operator.
    DPCell_ &
    operator=(DPCell_ const & other)
    {
        if (this != &other)
            _score = other._score;
            _horizontalScore = other._horizontalScore;
            _verticalScore = other._verticalScore;
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

// Needed for banded chain alignment for the set.
template <typename TScoreValueLeft, typename TScoreValueRight>
inline bool operator<(DPCell_<TScoreValueLeft, AffineGaps> const & left,
                      DPCell_<TScoreValueRight, AffineGaps> const & right)
{
    return left._score < right._score && left._horizontalScore < right._horizontalScore &&
           left._verticalScore < right._verticalScore;
}

// ----------------------------------------------------------------------------
// Function _verticalScoreOfCell()
// ----------------------------------------------------------------------------

// Returns the score of the matrix for vertical-gaps of the given cell.
template <typename TScoreValue>
inline typename  Reference<DPCell_<TScoreValue, AffineGaps> >::Type
_verticalScoreOfCell(DPCell_<TScoreValue, AffineGaps> & dpCell)
{
    return dpCell._verticalScore;
}

template <typename TScoreValue>
inline typename  Reference<DPCell_<TScoreValue, AffineGaps> const>::Type
_verticalScoreOfCell(DPCell_<TScoreValue, AffineGaps> const & dpCell)
{
    return dpCell._verticalScore;
}

// ----------------------------------------------------------------------------
// Function _setVerticalScoreOfCell()
// ----------------------------------------------------------------------------

// Returns the score of the matrix for vertical-gaps of the given cell.
template <typename TScoreValue>
inline void
_setVerticalScoreOfCell(DPCell_<TScoreValue, AffineGaps> & dpCell, TScoreValue const & newVerticalScore)
{
    dpCell._verticalScore = newVerticalScore;
}

// ----------------------------------------------------------------------------
// Function _horizontalScoreOfCell()
// ----------------------------------------------------------------------------

// Returns the score of the matrix for horizontal-gaps of the given cell.
template <typename TScoreValue>
inline typename  Reference<DPCell_<TScoreValue, AffineGaps> >::Type
_horizontalScoreOfCell(DPCell_<TScoreValue, AffineGaps> & dpCell)
{
    return dpCell._horizontalScore;
}

template <typename TScoreValue>
inline typename  Reference<DPCell_<TScoreValue, AffineGaps> const>::Type
_horizontalScoreOfCell(DPCell_<TScoreValue, AffineGaps> const & dpCell)
{
    return dpCell._horizontalScore;
}

// ----------------------------------------------------------------------------
// Function _setHorizontalScoreOfCell()
// ----------------------------------------------------------------------------

// Returns the score of the matrix for vertical-gaps of the given cell.
template <typename TScoreValue>
inline void
_setHorizontalScoreOfCell(DPCell_<TScoreValue, AffineGaps> & dpCell, TScoreValue const & newHorizontalScore)
{
    dpCell._horizontalScore = newHorizontalScore;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_CELL_AFFINE_H_
