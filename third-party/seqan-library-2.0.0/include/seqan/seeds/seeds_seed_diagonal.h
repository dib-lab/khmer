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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// SeedDiagonal objects store information about seed parts/segments of a
// ChainedSeed.
// ==========================================================================

#ifndef SEQAN_SEEDS_SEEDS_SEED_DIAGONAL_H_
#define SEQAN_SEEDS_SEEDS_SEED_DIAGONAL_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// ---------------------------------------------------------------------------
// Class SeedDiagonal
// ---------------------------------------------------------------------------

/*!
 * @class SeedDiagonal
 * @headerfile <seqan/seeds.h>
 * @brief Store the information about a seed segment.
 *
 * @signature template <typename TPosition, typename TSize>
 *            class SeedDiagonal;
 *
 * @tparam TSize The type to use for the seed length.
 * @tparam TPosition The type to use for positions.
 */

/*!
 * @fn SeedDiagonal::SeedDiagonal
 * @brief Constructor
 *
 * @signature SeedDiagonal::SeedDiagonal();
 * @signature SeedDiagonal::SeedDiagonal(beginPosH, beginPosV, length);
 *
 * @param[in] beginPosH The begin position in the horizontal orientation.
 * @param[in] beginPosV The begin position in the vertical orientation.
 * @param[in] length    The length of the seed diagonal.
 */

/*!
 * @var TPosition SeedDiagonal::beginPositionH
 * @brief The position in the database sequence (horizontal).
 *
 * @var TPosition SeedDiagonal::beginPositionV
 * @brief The position in the query sequence (vertical).
 *
 * @var TSize SeedDiagonal::length
 * @brief The length of the diagonal.
 */

template <typename TPosition, typename TSize>
class SeedDiagonal
{
public:
    TPosition beginPositionH;
    TPosition beginPositionV;
    TSize length;

    SeedDiagonal() : beginPositionH(0), beginPositionV(0), length(0)
    {}

    SeedDiagonal(TPosition _beginPositionH, TPosition _beginPositionV, TSize _length) :
            beginPositionH(_beginPositionH), beginPositionV(_beginPositionV), length(_length)
    {}
};

// ===========================================================================
// Metafunctions
// ===========================================================================

// ---------------------------------------------------------------------------
// Metafunction Position
// ---------------------------------------------------------------------------

/*!
 * @mfn SeedDiagonal#Position
 * @brief The position type of a @link SeedDiagonal @endlink.
 *
 * @signature Position<TDiagonal>::Type;
 *
 * @tparam TDiagonal A @link SeedDiagonal @endlink.
 *
 * @return Type The position type of <tt>TDiagonal</tt>.
 */

template <typename TPosition, typename TSize>
struct Position<SeedDiagonal<TPosition, TSize> >
{
    typedef TPosition Type;
};

template <typename TPosition, typename TSize>
struct Position<SeedDiagonal<TPosition, TSize> const> : Position<SeedDiagonal<TPosition, TSize> >
{};

// ---------------------------------------------------------------------------
// Metafunction Size
// ---------------------------------------------------------------------------

/*!
 * @mfn SeedDiagonal#Size
 * @brief The size type of a @link SeedDiagonal @endlink.
 *
 * @signature Size<TDiagonal>::Type;
 *
 * @tparam TDiagonal A @link SeedDiagonal @endlink.
 *
 * @return Type The size type of <tt>TDiagonal</tt>.
 */

template <typename TPosition, typename TSize>
struct Size<SeedDiagonal<TPosition, TSize> >
{
  typedef TSize Type;
};

template <typename TPosition, typename TSize>
struct Size<SeedDiagonal<TPosition, TSize> const> : Size<SeedDiagonal<TPosition, TSize> >
{};

// ===========================================================================
// Functions
// ===========================================================================

// ---------------------------------------------------------------------------
// Debugger Function operator<<()
// ---------------------------------------------------------------------------

template <typename TStream, typename TPosition, typename TSize>
inline TStream &
operator<<(TStream & stream, SeedDiagonal<TPosition, TSize> const & seedDiagonal)
{
    return stream << "SeedDiagonal(" << seedDiagonal.beginPositionH << ", " << seedDiagonal.beginPositionV
                  << ", "  << seedDiagonal.length << ")";
}

// ---------------------------------------------------------------------------
// Function operator==()
// ---------------------------------------------------------------------------

template <typename TPosition, typename TSize>
inline bool
operator==(SeedDiagonal<TPosition, TSize> const & a, SeedDiagonal<TPosition, TSize> const & b)
{
    return a.beginPositionH == b.beginPositionH &&
            a.beginPositionV == b.beginPositionV &&
            a.length == b.length;
}

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_SEED_DIAGONAL_H_
