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
// Declares the DPCell, which is used to substitute the score of each cell.
// Thus, we are able to add additional features to an alignment cell such as
// a flag to indicate if it is forbidden or not. Or to store two additional
// scores necessary for the affine gap function.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_CELL_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_CELL_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class DPCell_
// ----------------------------------------------------------------------------

// Used to store the score of a particular cell of the score matrix.
// It can be specialized for linear and affine gap costs.
// For affine gap costs it stores the values of all three matrices at a particular
// position of the matrix.
template <typename TScoreValue, typename TGapCosts>
class DPCell_;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TGapCostFunction>
struct Value<DPCell_<TScoreValue, TGapCostFunction> >
{
    typedef TScoreValue Type;
};

template <typename TScoreValue, typename TGapCostFunction>
struct Value<DPCell_<TScoreValue, TGapCostFunction> const>
{
    typedef TScoreValue const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TGapCostFunction>
struct Reference<DPCell_<TScoreValue, TGapCostFunction> >
{
    typedef TScoreValue & Type;
};

template <typename TScoreValue, typename TGapCostFunction>
struct Reference<DPCell_<TScoreValue, TGapCostFunction> const>
{
    typedef TScoreValue const & Type;
};

// ----------------------------------------------------------------------------
// Metafunction DPCellDefaultInfinity
// ----------------------------------------------------------------------------

// Defines the default infinity value for a DPCell.
template <typename T>
struct DPCellDefaultInfinity
{
    static const int VALUE;
};

template <typename T>
const int DPCellDefaultInfinity<T>::VALUE = MinValue<int>::VALUE;

// We use the min value of the score type and shift it one bits to the left.  This way we can use "infinity" without
// checking for it during the computation.
template <typename TScoreValue, typename TGapCostFunction>
struct DPCellDefaultInfinity<DPCell_<TScoreValue, TGapCostFunction> >
{
    static const TScoreValue VALUE;
};

template <typename TScoreValue, typename TGapCostFunction>
const TScoreValue DPCellDefaultInfinity<DPCell_<TScoreValue, TGapCostFunction> >::VALUE = MinValue<TScoreValue>::VALUE / 2;

template <typename TScoreValue, typename TGapCostFunction>
struct DPCellDefaultInfinity<DPCell_<TScoreValue, TGapCostFunction> const>:
    DPCellDefaultInfinity<DPCell_<TScoreValue, TGapCostFunction> >{};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _scoreOfCell
// ----------------------------------------------------------------------------

// Returns the score value for a given cell.
template <typename TScoreValue, typename TGapCosts>
inline typename Reference<DPCell_<TScoreValue, TGapCosts> >::Type
_scoreOfCell(DPCell_<TScoreValue, TGapCosts> & dpCell)
{
    return dpCell._score;
}

template <typename TScoreValue, typename TGapCosts>
inline typename Reference<DPCell_<TScoreValue, TGapCosts> const>::Type
_scoreOfCell(DPCell_<TScoreValue, TGapCosts> const & dpCell)
{
    return dpCell._score;
}

// ----------------------------------------------------------------------------
// Function _setScoreOfCell
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TGapCosts>
inline void
_setScoreOfCell(DPCell_<TScoreValue, TGapCosts> & dpCell, TScoreValue const & newScore)
{
    dpCell._score = newScore;
}

// ----------------------------------------------------------------------------
// Function _verticalScoreOfCell()
// ----------------------------------------------------------------------------

// Returns the score of the matrix for vertical-gaps of the given cell.
template <typename TScoreValue, typename TGapSpec>
inline typename  Reference<DPCell_<TScoreValue, TGapSpec> >::Type
_verticalScoreOfCell(DPCell_<TScoreValue, TGapSpec> & dpCell)
{
    return dpCell._score;
}

template <typename TScoreValue, typename TGapSpec>
inline typename  Reference<DPCell_<TScoreValue, TGapSpec> const>::Type
_verticalScoreOfCell(DPCell_<TScoreValue, TGapSpec> const & dpCell)
{
    return dpCell._score;
}

// ----------------------------------------------------------------------------
// Function _setVerticalScoreOfCell()
// ----------------------------------------------------------------------------

// Returns the score of the matrix for vertical-gaps of the given cell.
template <typename TScoreValue, typename TGapSpec>
inline void
_setVerticalScoreOfCell(DPCell_<TScoreValue, TGapSpec> & /*dpCell*/, TScoreValue const & /*newVerticalScore*/)
{
    // no-op
}

// ----------------------------------------------------------------------------
// Function _horizontalScoreOfCell()
// ----------------------------------------------------------------------------

// Returns the score of the matrix for horizontal-gaps of the given cell.
template <typename TScoreValue, typename TGapSpec>
inline typename  Reference<DPCell_<TScoreValue, TGapSpec> >::Type
_horizontalScoreOfCell(DPCell_<TScoreValue, TGapSpec> & dpCell)
{
    return dpCell._score;
}

template <typename TScoreValue, typename TGapSpec>
inline typename  Reference<DPCell_<TScoreValue, TGapSpec> const>::Type
_horizontalScoreOfCell(DPCell_<TScoreValue, TGapSpec> const & dpCell)
{
	return dpCell._score;
}

// ----------------------------------------------------------------------------
// Function _setHorizontalScoreOfCell()
// ----------------------------------------------------------------------------

// Returns the score of the matrix for vertical-gaps of the given cell.
template <typename TScoreValue, typename TGapSpec>
inline void
_setHorizontalScoreOfCell(DPCell_<TScoreValue, TGapSpec> & /*dpCell*/, TScoreValue const & /*newHorizontalScore*/)
{
    // no-op
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_CELL_H_
