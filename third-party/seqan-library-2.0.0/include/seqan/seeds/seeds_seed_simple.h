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
// Specialization "Simple" for class Seed.
// ==========================================================================

#ifndef SEQAN_SEEDS_SEEDS_SEED_SIMPLE_H_
#define SEQAN_SEEDS_SEEDS_SEED_SIMPLE_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// ---------------------------------------------------------------------------
// Class Simple Seed
// ---------------------------------------------------------------------------

/*!
 * @class SimpleSeed
 * @headerfile <seqan/seeds.h>
 * @extends Seed
 * @brief Describes a seed with start and end position and diagonal upper and lower bounds.
 *
 * @signature template <[typename TConfig]>
 *            class Seed<Simple, TConfig>;
 *
 * @tparam TConfig The configuration to use.  Defaults to @link DefaultSeedConfig @endlink.
 */

/*!
 * @fn SimpleSeed::Seed
 * @brief Constructor
 *
 * @signature Seed::Seed();
 * @signature Seed::Seed(beginPosH, beginPosV, length);
 * @signature Seed::Seed(beginPosH, beginPosV, endPosH, endPosV);
 *
 * @param[in] beginPosH The begin position in the horizontal position.
 * @param[in] beginPosV The begin position in the vertical position.
 * @param[in] length    The length of the seed (in both directions).
 * @param[in] endPosH   The end position in the horizontal position.
 * @param[in] endPosV   The end position in the vertical position.
 */

template <typename TConfiguration>
class Seed<Simple, TConfiguration>
{
public:
    typedef typename TConfiguration::TPosition TPosition;
    typedef typename TConfiguration::TSize TSize;
    typedef typename TConfiguration::TDiagonal TDiagonal;
    typedef typename TConfiguration::TScoreValue TScoreValue;

    TPosition _beginPositionH;
    TPosition _beginPositionV;
    TPosition _endPositionH;
    TPosition _endPositionV;
    TDiagonal _lowerDiagonal;
    TDiagonal _upperDiagonal;
    TScoreValue _score;

    Seed() :_beginPositionH(0), _beginPositionV(0), _endPositionH(0), _endPositionV(0),_lowerDiagonal(0), _upperDiagonal(0),
            _score()
    {}

    Seed(TPosition beginPositionH, TPosition beginPositionV, TPosition seedLength) :
            _beginPositionH(beginPositionH), _beginPositionV(beginPositionV), _endPositionH(beginPositionH + seedLength),
            _endPositionV(beginPositionV + seedLength), _lowerDiagonal(static_cast<TDiagonal>(beginPositionH - beginPositionV)),
            _upperDiagonal(static_cast<TDiagonal>(beginPositionH - beginPositionV)), _score(0)
    {
        SEQAN_ASSERT_GEQ(_upperDiagonal, _lowerDiagonal);
    }

    Seed(TPosition beginPositionH, TPosition beginPositionV, TPosition endPositionH, TPosition endPositionV) :
            _beginPositionH(beginPositionH),
            _beginPositionV(beginPositionV),
            _endPositionH(endPositionH),
            _endPositionV(endPositionV),
            _lowerDiagonal(std::min(static_cast<TDiagonal>(beginPositionH - beginPositionV),
                                    static_cast<TDiagonal>(endPositionH - endPositionV))),
            _upperDiagonal(std::max(static_cast<TDiagonal>(beginPositionH - beginPositionV),
                                    static_cast<TDiagonal>(endPositionH - endPositionV))),
            _score(0)
    {
        SEQAN_ASSERT_GEQ(_upperDiagonal, _lowerDiagonal);
    }

    template <typename TSeed2>
    Seed(TSeed2 const & other) :
              _beginPositionH(beginPositionH(other)),
              _beginPositionV(beginPositionV(other)),
              _endPositionH(endPositionH(other)),
              _endPositionV(endPositionV(other)),
              _lowerDiagonal(lowerDiagonal(other)),
              _upperDiagonal(upperDiagonal(other)),
              _score(0)
    {
        SEQAN_ASSERT_GEQ(_upperDiagonal, _lowerDiagonal);
    }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

// ---------------------------------------------------------------------------
// Debug Function operator<<()
// ---------------------------------------------------------------------------

template <typename TStream, typename TConfig>
inline TStream &
operator<<(TStream & stream, Seed<Simple, TConfig> const & seed)
{
    return stream << "Seed<Simple, TConfig>(" << beginPositionH(seed)
                  << ", " << beginPositionV(seed) << ", "
                  << endPositionH(seed) << ", "
                  << endPositionV(seed) << ", lower diag = "
                  << lowerDiagonal(seed) << ", upper diag = "
                  << upperDiagonal(seed) << ")";
}

// ---------------------------------------------------------------------------
// Function operator==()
// ---------------------------------------------------------------------------

template <typename TConfig>
inline bool
operator==(Seed<Simple, TConfig> const & a, Seed<Simple, TConfig> const & b)
{
    return a._beginPositionH == b._beginPositionH &&
            a._beginPositionV == b._beginPositionV &&
            a._endPositionH == b._endPositionH &&
            a._endPositionV == b._endPositionV &&
            a._upperDiagonal == b._upperDiagonal &&
            a._lowerDiagonal == b._lowerDiagonal;
}

// ---------------------------------------------------------------------------
// Function beginPositionH()
// ---------------------------------------------------------------------------

template <typename TConfig>
inline typename Position<Seed<Simple, TConfig> >::Type
beginPositionH(Seed<Simple, TConfig> const & seed)
{
    return seed._beginPositionH;
}

// ---------------------------------------------------------------------------
// Function setBeginPositionH()
// ---------------------------------------------------------------------------

template <typename TConfig, typename TPosition>
inline void
setBeginPositionH(Seed<Simple, TConfig> & seed, TPosition pos)
{
    seed._beginPositionH = pos;
}

// ---------------------------------------------------------------------------
// Function endPositionH()
// ---------------------------------------------------------------------------

template <typename TConfig>
inline typename Position<Seed<Simple, TConfig> >::Type
endPositionH(Seed<Simple, TConfig> const & seed)
{
    return seed._endPositionH;
}

// ---------------------------------------------------------------------------
// Function setEndPositionH()
// ---------------------------------------------------------------------------

template <typename TConfig, typename TPosition>
inline void
setEndPositionH(Seed<Simple, TConfig> & seed, TPosition pos)
{
    seed._endPositionH = pos;
}

// ---------------------------------------------------------------------------
// Function beginPositionV()
// ---------------------------------------------------------------------------

template <typename TConfig>
inline typename Position<Seed<Simple, TConfig> >::Type
beginPositionV(Seed<Simple, TConfig> const & seed)
{
    return seed._beginPositionV;
}

// ---------------------------------------------------------------------------
// Function setBeginPositionV()
// ---------------------------------------------------------------------------

template <typename TConfig, typename TPosition>
inline void
setBeginPositionV(Seed<Simple, TConfig> & seed, TPosition pos)
{
    seed._beginPositionV = pos;
}

// ---------------------------------------------------------------------------
// Function endPositionV()
// ---------------------------------------------------------------------------

template <typename TConfig>
inline typename Position<Seed<Simple, TConfig> >::Type
endPositionV(Seed<Simple, TConfig> const & seed)
{
    return seed._endPositionV;
}

// ---------------------------------------------------------------------------
// Function setEndPositionV()
// ---------------------------------------------------------------------------

template <typename TConfig, typename TPosition>
inline void
setEndPositionV(Seed<Simple, TConfig> & seed, TPosition pos)
{
    seed._endPositionV = pos;
}

// ---------------------------------------------------------------------------
// Debug Function __write()
// ---------------------------------------------------------------------------

struct Tikz_ {};

template <typename TStream, typename TConfig>
inline void
__write(TStream & stream, Seed<Simple, TConfig> const & seed, Tikz_ const &)
{
//IOREV _nodoc_ specialization not documented
    stream << "\\draw[seed] (" << getBeginDim1(seed) << ", -" << getBeginDim0(seed) << ") -- (" << (getEndDim1(seed) - 1) << ", -" << (getEndDim0(seed) - 1) << ");" << std::endl;
}

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_SEED_SIMPLE_H_
