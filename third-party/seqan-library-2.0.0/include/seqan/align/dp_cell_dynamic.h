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
// Implements the dynamic gap model published in "Dynamic Gaps Selector:
// A Smith Waterman Sequence Alignment Algorithm with Affine Gap Model
// Optimization" by Gianvito Urgese et al.
// ==========================================================================

#ifndef INCLUDE_SEQAN_ALIGN_DP_CELL_DYNAMIC_H_
#define INCLUDE_SEQAN_ALIGN_DP_CELL_DYNAMIC_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct DynamicGapExtensionHorizontal_;
typedef Tag<DynamicGapExtensionHorizontal_> DynamicGapExtensionHorizontal;

struct DynamicGapExtensionVertical_;
typedef Tag<DynamicGapExtensionVertical_> DynamicGapExtensionVertical;

enum DynamicGapsMask
{
    MASK_VERTICAL_GAP = 1,
    MASK_HORIZONTAL_GAP = 2
};

// ----------------------------------------------------------------------------
// Class DPCell                                                   [DynamicGaps]
// ----------------------------------------------------------------------------

// The specialization for linear gap cost function.
// It solely stores the maximal score.
template <typename TScoreValue>
class DPCell_<TScoreValue, DynamicGaps>
{
public:
    TScoreValue _score;
    char _flagMask;

    // The default c'tor.
    DPCell_() : _score(DPCellDefaultInfinity<DPCell_>::VALUE), _flagMask(0)
    {}

    // The copy c'tor.
    DPCell_(DPCell_<TScoreValue, DynamicGaps> const & other) : _score(other._score), _flagMask(other._flagMask)
    {}

    // Implicit c'tor.
    DPCell_(TScoreValue const & score) : _score(score), _flagMask(0)
    {}

    // The assignment operator.
    DPCell_ &
    operator=(DPCell_<TScoreValue, DynamicGaps> const & other)
    {
        if (this != &other)
        {
            _score = other._score;
            _flagMask = other._flagMask;
        }
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

template <typename TCell, typename TBoolV, typename TBoolH>
struct SetGapExtension;

template <typename TScoreValue>
struct SetGapExtension<DPCell_<TScoreValue, DynamicGaps>, False, False>
{
    static const char VALUE = 0;
};

template <typename TScoreValue>
struct SetGapExtension<DPCell_<TScoreValue, DynamicGaps>, False, True>
{
    static const char VALUE = MASK_HORIZONTAL_GAP;
};

template <typename TScoreValue>
struct SetGapExtension<DPCell_<TScoreValue, DynamicGaps>, True, False>
{
    static const char VALUE = MASK_VERTICAL_GAP;
};

template <typename TScoreValue>
struct SetGapExtension<DPCell_<TScoreValue, DynamicGaps>, True, True>
{
    static const char VALUE = MASK_HORIZONTAL_GAP | MASK_VERTICAL_GAP;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TScoreValue, typename TFlag>
inline void _setBit(DPCell_<TScoreValue, DynamicGaps> & cell,
                    TFlag const & /*flag*/,
                    DynamicGapExtensionVertical const & /*tag*/)
{
    if (IsSameType<TFlag, True>::VALUE)
        cell._flagMask |= MASK_VERTICAL_GAP;
    else
        cell._flagMask &= ~MASK_VERTICAL_GAP;
}

template <typename TScoreValue, typename TFlag>
inline void _setBit(DPCell_<TScoreValue, DynamicGaps> & cell,
                    TFlag const & /*flag*/,
                    DynamicGapExtensionHorizontal const & /*tag*/)
{
    if (IsSameType<TFlag, True>::VALUE)
        cell._flagMask |= MASK_HORIZONTAL_GAP;
    else
        cell._flagMask &= ~MASK_HORIZONTAL_GAP;
}

template <typename TScoreValue, typename TSpec>
inline bool isGapExtension(DPCell_<TScoreValue, DynamicGaps> const & cell,
                           TSpec const & /*spec*/)
{
    if (IsSameType<TSpec, DynamicGapExtensionHorizontal>::VALUE)
        return cell._flagMask & MASK_HORIZONTAL_GAP;
    else
        return cell._flagMask & MASK_VERTICAL_GAP;
}

template <typename TScoreValue, typename TFlagV, typename TFlagH>
inline void
setGapExtension(DPCell_<TScoreValue, DynamicGaps> & cell,
                TFlagV const & /*vert*/,
                TFlagH const & /*hori*/)
{
    cell._flagMask = SetGapExtension<DPCell_<TScoreValue, DynamicGaps>, TFlagV, TFlagH>::VALUE;
}

// ----------------------------------------------------------------------------
// Function operator<()
// ----------------------------------------------------------------------------

// Needed for banded chain alignment for the std::set.
template <typename TScoreValueLeft, typename TScoreValueRight>
inline bool operator<(DPCell_<TScoreValueLeft, DynamicGaps> const & left,
                      DPCell_<TScoreValueRight, DynamicGaps> const & right)
{
    return left._score < right._score;
}

}  // namespace seqan

#endif  // INCLUDE_SEQAN_ALIGN_DP_CELL_DYNAMIC_H_
