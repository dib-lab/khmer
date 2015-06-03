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
// This class stores the band information as well as the meta-inforation,
// whether a band was selected or not.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_BAND_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_DP_BAND_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @defgroup DPBandSwitch
 * @brief Tags used to switch between banded and unbanded alignment.
 *
 * @tag DPBandSwitch#BandOn
 * @brief Switches banded alignment on.
 * @headerfile <seqan/align.h>
 * @signature struct BandOn_;
 * @signature typedef Tag<BandOn_> BandOn;
 *
 * @tag DPBandSwitch#BandOff
 * @brief Switches banded alignment off.
 * @headerfile <seqan/align.h>
 * @signature struct BandOff_;
 * @signature typedef Tag<BandOff_> BandOff;
 */

// ----------------------------------------------------------------------------
// Tag BandOff
// ----------------------------------------------------------------------------

// Used when computing unbanded alignments.
struct BandOff_;
typedef Tag<BandOff_> BandOff;

// ----------------------------------------------------------------------------
// Tag BandOn
// ----------------------------------------------------------------------------

// Used when computing banded alignments.
struct BandOn_;
typedef Tag<BandOn_> BandOn;

// ----------------------------------------------------------------------------
// Class DPBandConfig
// ----------------------------------------------------------------------------

/*!
 * @class DPBandConfig
 * @headerfile <seqan/align.h>
 * @brief Simple class to configure banded alignments.
 *
 * @signature template <typename TSwitch>
 *            class DPBandConfig<TSwitch>;
 *
 * @tparam TSwitch Tag to switch between banded and unbanded alignments.
 *              One of @link DPBandSwitch @endlink. Defaults to @link DPBandSwitch#BandOff @endlink.
 *
 * To compute banded alignments use @link DPBand @endlink as a shortcut for the <tt>DPBandConfig</tt> with
 * band switched on.
 */

// Simple band class.
template <typename TSpec = BandOff>
struct DPBandConfig {};

// ----------------------------------------------------------------------------
// Class DPBandConfig                                                  [BandOff]
// ----------------------------------------------------------------------------

// The specialization when using no band.
// Per default the member variables _lowerDiagonal and _upperDiagonal are
// always 0.
template <>
struct DPBandConfig<BandOff>
{
    typedef int TPosition;
};

// ----------------------------------------------------------------------------
// Class DPBandConfig                                                  [BandOn]
// ----------------------------------------------------------------------------

// The specialization when using a band.
// On construction the diagonals are set to 0.
template <>
struct DPBandConfig<BandOn>
{
    typedef int TPosition;

    int _lowerDiagonal;
    int _upperDiagonal;

/*!
 * @fn DPBandConfig::DPBandConfig
 * @brief Constructor.
 *
 * @signature DPBandConfig<TSwitch>();
 * @signature DPBandConfig<BandOn>(lowerDiag, upperDiag);
 *
 * @tparam TSwitch Tag to switch between banded and unbanded alignments. One of @link DPBandSwitch @endlink.
 *              The second constructor is only supported when @link DPBandConfig @endlink is specialized with
 *              @link DPBandSwitch#BandOn @endlink.
 *
 * @param lowerDiag The value for the lower diagonal of the band.
 * @param upperDiag The value for the upper diagonal of the band.
 *
 * A negative value for the diagonals indicates an intersection of the diagonal with the vertical sequence (y-axis)
 * and a positive value indicates an intersection with the horizontal sequence (x-axis).
 * The value of the lower diagonal has to compare less or equal to the value of the upper diagonal.
 */

    DPBandConfig() :
        _lowerDiagonal(0), _upperDiagonal(0) {}

    DPBandConfig(int lowerDiagonal, int upperDiagonal) :
        _lowerDiagonal(lowerDiagonal), _upperDiagonal(upperDiagonal)
    {
        SEQAN_ASSERT_LEQ(lowerDiagonal, upperDiagonal);
    }
};

/*!
 * @typedef DPBand
 * @headerfile <seqan/align.h>
 * @brief Global typedef used for @link DPBandConfig @endlink specialized with @link DPBandSwitch#BandOn @endlink.
 *
 * @signature typedef DPBandConfig<BandOn> DPBand;
 */

// Typedef for Band.
typedef DPBandConfig<BandOn> DPBand;

// ============================================================================
// Metafunctions
// ============================================================================

/*!
 * @mfn DPBandConfig#Position
 * @headerfile <seqan/align.h>
 * @brief Metafunction returning the position type.
 *
 * @signature typename Position<T>::Type;
 *
 * @tparam T The type @link DPBandConfig @endlink to query the position type for.
 * @return TPosition The position type.
 */

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename T>
struct Position<DPBandConfig<T> >
{
    typedef typename DPBandConfig<T>::TPosition Type;
};

template <typename T>
struct Position<DPBandConfig<T> const>:
    Position<DPBandConfig<T> >{};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

/*!
 * @mfn DPBandConfig#Size
 * @headerfile <seqan/align.h>
 * @brief Metafunction returning the size type.
 *
 * @signature typename Size<T>::Type;
 *
 * @tparam T The type @link DPBandConfig @endlink to query the size type for.
 * @return TSize The size type.
 */

template <typename T>
struct Size<DPBandConfig<T> >
{
    typedef unsigned Type;
};


template <typename T>
struct Size<DPBandConfig<T> const>:
    Size<DPBandConfig<T> >{};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setLowerDiagonal
// ----------------------------------------------------------------------------

/*!
 * @fn DPBandConfig#setLowerDiagonal
 * @headerfile <seqan/align.h>
 * @brief Sets the value of the lower diagonal.
 *
 * @signature setLowerDiagonal(obj, val);
 *
 * @param obj The object of type @link DPBandConfig @endlink to set the lower diagonal for.
 * @param val The new value for the lower diagonal.
 *
 * @note If the band is switched off, this function defaults to no-op.
 */

inline void
setLowerDiagonal(DPBandConfig<BandOff> & /*dpBand*/, int /*newLowerDiagonal*/)
{
    //no-op
}

inline void
setLowerDiagonal(DPBandConfig<BandOn> & dpBand, int newLowerDiagonal)
{
    dpBand._lowerDiagonal = newLowerDiagonal;
}

// ----------------------------------------------------------------------------
// Function lowerDiagonal
// ----------------------------------------------------------------------------

/*!
 * @fn DPBandConfig#lowerDiagonal
 * @headerfile <seqan/align.h>
 * @brief Returns the value of the lower diagonal.
 *
 * @signature TPosition lowerDiagonal(obj);
 *
 * @param obj The object of type @link DPBandConfig @endlink to query the lower diagonal for.
 *
 * @note If the band is switched off this function always returns 0.
 * @return TPosition The value of the lower diagonal.
 */

inline Position<DPBandConfig<BandOff> >::Type
lowerDiagonal(DPBandConfig<BandOff> const & /*dpBand*/)
{
    return 0;
}

template <typename TSwitch>
inline typename Position<DPBandConfig<BandOff> >::Type
lowerDiagonal(DPBandConfig<TSwitch> const & dpBand)
{
    return dpBand._lowerDiagonal;
}

// ----------------------------------------------------------------------------
// Function setUpperDiagonal
// ----------------------------------------------------------------------------

/*!
 * @fn DPBandConfig#setUpperDiagonal
 * @headerfile <seqan/align.h>
 * @brief Sets the value of the upper diagonal.
 *
 * @signature setUpperDiagonal(obj, val);
 *
 * @param obj The object of type @link DPBandConfig @endlink to set the upper diagonal for.
 * @param val The new value for the upper diagonal.
 *
 * @note If the band is switched off, this function defaults to no-op.
 */

inline void
setUpperDiagonal(DPBandConfig<BandOff> & /*dpBand*/, int /*newUpperDiagonal*/)
{
    //no-op
}

inline void
setUpperDiagonal(DPBandConfig<BandOn> & dpBand, int newUpperDiagonal)
{
    dpBand._upperDiagonal = newUpperDiagonal;
}

// ----------------------------------------------------------------------------
// Function upperDiagonal
// ----------------------------------------------------------------------------

/*!
 * @fn DPBandConfig#upperDiagonal
 * @headerfile <seqan/align.h>
 * @brief Returns the value of the upper diagonal.
 *
 * @signature TPosition upperDiagonal(obj);
 *
 * @param obj The object of type @link DPBandConfig @endlink to query the upper diagonal for.
 *
 * @note If the band is switched off this function always returns 0.
 * @return TPosition The value of the upper diagonal.
 */

inline Position<DPBandConfig<BandOff> >::Type
upperDiagonal(DPBandConfig<BandOff> const & /*dpBand*/)
{
    return 0;
}

template <typename TSwitch>
inline typename Position<DPBandConfig<BandOn> >::Type
upperDiagonal(DPBandConfig<TSwitch> const & dpBand)
{
    return dpBand._upperDiagonal;
}

// ----------------------------------------------------------------------------
// Function bandSize()
// ----------------------------------------------------------------------------

/*!
 * @fn DPBandConfig#bandSize
 * @headerfile <seqan/align.h>
 * @brief Returns the size of the band.
 *
 * @signature TSize bandSize(obj);
 *
 * @param obj The object of type @link DPBandConfig @endlink to query the band size for.
 *
 * @note If the band is switched off this function always returns 0.
 * @return TSize The number of diagonals covered by the band.
 */

inline Size<DPBandConfig<BandOff> >::Type
bandSize(DPBandConfig<BandOff> const &)
{
    return 0;
}

template <typename TSwitch>
inline typename Size<DPBandConfig<TSwitch> >::Type
bandSize(DPBandConfig<TSwitch> const & band)
{
    return upperDiagonal(band) - lowerDiagonal(band) + 1;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_BAND_H_
