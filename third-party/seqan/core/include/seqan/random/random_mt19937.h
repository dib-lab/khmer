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

/**
.Spec.Mersenne Twister Rng
..general:Class.Rng
..signature:Rng<MersenneTwister>
..summary:Mersenne Twister 19937 Random Number Generator
..cat:Random
..include:seqan/random.h
..wiki:Tutorial/Randomness|Tutorial: Randomness
*/
template <>
class Rng<MersenneTwister>
{
public:
    ext::MTRand _mtRand;


/**
.Memfunc.Mersenne Twister Rng#Rng
..class:Spec.Mersenne Twister Rng
..summary:Constructor Mersenne Twister Rng.
..signature:Rng<MersenneTwister>([seed])
..param.seed:Seed for the initialization of the Mersenne Twister, defaults to 0.
...type:nolink:double.
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

/**
.Function.reSeed
..class:Spec.Mersenne Twister Rng
..summary:Reset and re-seed MersenneTwister
..cat:Random
..signature:reSeed(mt[, seed])
..param.mt:The @Spec.Mersenne Twister Rng@ to reset.
...type:Spec.Mersenne Twister Rng
..param.seed:Optional seed to initialize the RNG with.
...default:0
...type:nolink:$__uint32$
..include:seqan/random.h
..wiki:Tutorial/Randomness|Tutorial: Randomness
*/

inline void
reSeed(Rng<MersenneTwister> & mt, __uint32 const seed)
{
    SEQAN_CHECKPOINT;

    mt._mtRand.seed(seed);
}

inline void
reSeed(Rng<MersenneTwister> & mt)
{
    SEQAN_CHECKPOINT;

    reSeed(mt, 0);
}

}  // namespace seqan

#endif  // SEQAN_RANDOM_RANDOM_MT19937_H_
