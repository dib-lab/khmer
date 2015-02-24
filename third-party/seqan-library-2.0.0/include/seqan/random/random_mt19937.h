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
// An implementation of the Mersenne Twister 19937 random number generator.
// ==========================================================================

#ifndef SEQAN_RANDOM_RANDOM_MT19937_H_
#define SEQAN_RANDOM_RANDOM_MT19937_H_

#include <seqan/random/ext_MersenneTwister.h>

namespace seqan {

// ===========================================================================
// Forwards, Tags.
// ===========================================================================

// Tag for selecting a mersenne twister.
struct MersenneTwister {};

// ===========================================================================
// Classes
// ===========================================================================

/*!
 * @class MersenneTwisterRng
 * @extends Rng
 * @headerfile <seqan/random.h>
 * @brief Mersenne Twister 19937 Random Number Generator.
 *
 * @signature template <>
 *            class Rng<MersenneTwister>;
 */

template <>
class Rng<MersenneTwister>
{
public:
    ext::MTRand _mtRand;

/*!
 * @fn MersenneTwisterRng::Rng
 * @brief Constructor Mersenne Twister Rng.
 *
 * @signature Rng::Rng([seed]);
 *
 * @param[in] seed The <tt>unsigned</tt> value to use for seeding, defaults to 0.
 */

    Rng() : _mtRand(0lu)
    { SEQAN_CHECKPOINT; }

    Rng(unsigned seed) : _mtRand(seed)
    { SEQAN_CHECKPOINT; }

    inline
    unsigned
    operator()()
    {
        SEQAN_CHECKPOINT;
        return _mtRand.randInt();
    }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

template <>
struct Value<Rng<MersenneTwister> >
{
    typedef unsigned Type;
};

template <>
struct Value<const Rng<MersenneTwister> > : Value<Rng<MersenneTwister> > {};

// ===========================================================================
// Functions
// ===========================================================================

inline unsigned
pickRandomNumber(Rng<MersenneTwister> & mt)
{
    SEQAN_CHECKPOINT;

    return mt._mtRand.randInt();
}

/*!
 * @fn MersenneTwisterRng#reSeed
 * @brief Reset and re-seed Mersenne Twister Rng.
 *
 * @signature void reSeed(mt[, seed]);
 *
 * @param[in,out] mt   The MersenneTwisterRng to re-seed.
 * @param[in]     seed The <tt>unsigned</tt> to use for re-seeding, defaults to 0.
 */

inline void
reSeed(Rng<MersenneTwister> & mt, __uint32 const seed = 0)
{
    SEQAN_CHECKPOINT;

    mt._mtRand.seed(seed);
}

}  // namespace seqan

#endif  // SEQAN_RANDOM_RANDOM_MT19937_H_
