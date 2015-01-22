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
// This class stores the band information as well as the meta-inforation,
// whether a band was selected or not.
// ==========================================================================

// TODO(holtgrew): Documentation in this header necessary or internal only?

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_BAND_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_BAND_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

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
// Class DPBand_
// ----------------------------------------------------------------------------

// Simple band class.
template <typename TSpec>
struct DPBand_ {};

// ----------------------------------------------------------------------------
// Class DPBand_                                                      [BandOff]
// ----------------------------------------------------------------------------

// The specialization when using no band.
// Per default the member variables _lowerDiagonal and _upperDiagonal are
// always 0.
template <>
struct DPBand_<BandOff>
{
    typedef int TPosition;
};

// ----------------------------------------------------------------------------
// Class DPBand_                                                       [BandOn]
// ----------------------------------------------------------------------------

// The specialization when using a band.
// On construction the diagonals are set to 0.
template <>
struct DPBand_<BandOn>
{
    typedef int TPosition;

    int _lowerDiagonal;
    int _upperDiagonal;

    DPBand_() :
        _lowerDiagonal(0), _upperDiagonal(0) {}

    DPBand_(int lowerDiagonal, int upperDiagonal) :
        _lowerDiagonal(lowerDiagonal), _upperDiagonal(upperDiagonal) {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename T>
struct Position<DPBand_<T> >
{
    typedef typename DPBand_<T>::TPosition Type;
};

template <typename T>
struct Position<DPBand_<T> const>:
    Position<DPBand_<T> >{};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename T>
struct Size<DPBand_<T> >:
    Position<DPBand_<T> >{};

template <typename T>
struct Size<DPBand_<T> const>:
    Size<DPBand_<T> >{};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setLowerDiagonal
// ----------------------------------------------------------------------------

inline void
setLowerDiagonal(DPBand_<BandOff> & /*dpBand*/, int /*newLowerDiagonal*/)
{
    //no-op
}

inline void
setLowerDiagonal(DPBand_<BandOn> & dpBand, int newLowerDiagonal)
{
    dpBand._lowerDiagonal = newLowerDiagonal;
}

// ----------------------------------------------------------------------------
// Function lowerDiagonal
// ----------------------------------------------------------------------------

inline int
lowerDiagonal(DPBand_<BandOff> const & /*dpBand*/)
{
    return 0;
}

template <typename TBandSwitch>
inline int
lowerDiagonal(DPBand_<TBandSwitch> const & dpBand)
{
    return dpBand._lowerDiagonal;
}

// ----------------------------------------------------------------------------
// Function setUpperDiagonal
// ----------------------------------------------------------------------------

inline void
setUpperDiagonal(DPBand_<BandOff> & /*dpBand*/, int /*newUpperDiagonal*/)
{
    //no-op
}

inline void
setUpperDiagonal(DPBand_<BandOn> & dpBand, int newUpperDiagonal)
{
    dpBand._upperDiagonal = newUpperDiagonal;
}

// ----------------------------------------------------------------------------
// Function upperDiagonal
// ----------------------------------------------------------------------------

inline int
upperDiagonal(DPBand_<BandOff> const & /*dpBand*/)
{
    return 0;
}

inline int
upperDiagonal(DPBand_<BandOn> const & dpBand)
{
    return dpBand._upperDiagonal;
}

// ----------------------------------------------------------------------------
// Function bandSize()
// ----------------------------------------------------------------------------

inline unsigned int
bandSize(DPBand_<BandOff> const &)
{
    return 0;
}

inline unsigned int
bandSize(DPBand_<BandOn> const & band)
{
    return upperDiagonal(band) - lowerDiagonal(band) + 1;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_BAND_H_
