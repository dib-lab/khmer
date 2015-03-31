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
// Basic definitions for the module random.
// ==========================================================================

#ifndef SEQAN_RANDOM_RANDOM_BASE_H_
#define SEQAN_RANDOM_RANDOM_BASE_H_

namespace seqan {

// ===========================================================================
// Forwards, Tags.
// ===========================================================================

// Forward to MersenneTwister, really defined in random_mt19937.h.
struct MersenneTwister;

/**
.Class.Rng:
..summary:Random Number Generator
..signature:Rng<>
..signature:Rng<TSpec>
..cat:Random
..param.TSpec:Random Number Generator specialization.
...default:@Spec.Mersenne Twister Rng@
..include:seqan/random.h
..wiki:Tutorial/Randomness|Tutorial: Randomness
..example
...text:The following code shows how to generate random numbers and shuffle a text.
...file:demos/random/random.cpp
...output:pickRandomNumber(rng) == 1608637542
pickRandomNumber(rng, uniformDouble) == 0.950714
pickRandomNumber(rng, uniformInt) == 27
pickRandomNumber(rng, normal) == 0.419823
pickRandomNumber(rng, logNormal) == 1.22431
pickRandomNumber(rng, logNormal2) == 2.78004
pickRandomNumber(rng, logNormal3) == 0.00155248
shuffle("Hello World!") ==  o!reWlloHld
*/

template <typename TSpec = MersenneTwister>
class Rng;

/**
.Class.Pdf:
..summary:ProbabilityDensityFunction
..signature:Pdf<TSpec>
..cat:Random
..param.TSpec:Specialization.
..include:seqan/random.h
..wiki:Tutorial/Randomness|Tutorial: Randomness
*/

template <typename TSpec>
class Pdf;

// ===========================================================================
// Classes
// ===========================================================================

/**
.Memfunc.Rng#operator()
..class:Class.Rng
..summary:Function call operator.
..signature:operator()
*/

// ===========================================================================
// Metafunctions
// ===========================================================================

///.Metafunction.Value.param.T.type:Class.Pdf
///.Metafunction.Value.class:Class.Pdf
// specification only

///.Metafunction.Value.param.T.type:Class.Rng
///.Metafunction.Value.class:Class.Rng
///.Metafunction.MinValue.param.T.type:Class.Rng
///.Metafunction.MinValue.class:Class.Rng
///.Metafunction.MaxValue.param.T.type:Class.Rng
///.Metafunction.MaxValue.class:Class.Rng

template <typename TSpec>
struct MaxValue<Rng<TSpec> >
{
    typedef typename Value<Rng<TSpec> >::Type TValue_;
    static const TValue_ VALUE;
};

template <typename TSpec>
const typename Value<Rng<TSpec> >::Type MaxValue<Rng<TSpec> >::VALUE = MaxValue<typename Value<Rng<TSpec> >::Type>::VALUE;

template <typename TSpec>
struct MaxValue<Rng<TSpec> const>
{
	typedef typename Value<Rng<TSpec> const>::Type TValue_;
	static const TValue_ VALUE;
};

template <typename TSpec>
const typename Value<Rng<TSpec> const>::Type MaxValue<Rng<TSpec> const>::VALUE = MaxValue<typename Value<Rng<TSpec> const>::Type>::VALUE;

template <typename TSpec>
struct MinValue<Rng<TSpec> >
{
	typedef typename Value<Rng<TSpec> >::Type TValue_;
	static const TValue_ VALUE;
};

template <typename TSpec>
const typename Value<Rng<TSpec> >::Type MinValue<Rng<TSpec> >::VALUE = MinValue<typename Value<Rng<TSpec> >::Type>::VALUE;

template <typename TSpec>
struct MinValue<Rng<TSpec> const>
{
	typedef typename Value<Rng<TSpec> const>::Type TValue_;
	static const TValue_ VALUE;
};

template <typename TSpec>
const typename Value<Rng<TSpec> const>::Type MinValue<Rng<TSpec> const>::VALUE = MinValue<typename Value<Rng<TSpec> const>::Type>::VALUE;

/**
.Metafunction.GetDefaultRng
..cat:Random
..summary:Return the default @Class.Rng|Random Number Generator@ to use in a given class, spezialiation or algorithm.
..signature:GetDefaultRng<T>::Type
..param.T:The class or algorithm tag to get the default @Class.Rng@ for.
..returns:The type of Rng specialization to use.
..remarks:Currently, the default value is @Spec.Mersenne Twister Rng@.
..see:Function.defaultRng
 */
template <typename T>
struct GetDefaultRng
{
    typedef Rng<MersenneTwister> Type;
};

// ===========================================================================
// Functions
// ===========================================================================

/**
.Function.pickRandomNumber
..class:Class.Rng
..class:Class.Pdf
..summary:Pick a random number using a random number generator object, possibly following the given distribution.
..cat:Random
..include:seqan/random.h
..wiki:Tutorial/Randomness|Tutorial: Randomness
..signature:pickRandomNumber(rng[, pdf])
..param.rng:Random number generator to use.
...type:Class.Rng
..param.pdf:Probability density function to use, if any.
...type:Class.Pdf
..returns:Random number as specified in pdf, if any, or rng. For more details refer to the SeqAn Tutorial.
 */
// specification only

/**
.Function.defaultRng
..summary:Default default random number generator object of a given type.
..cat:Random
..signature:defaultRng<TRng>()
..param.TRng:Type of the @Class.Rng@ to return the global default object of.
...default:$MersenneTwister$
...type:nolink:$MersenneTwister$
..returns:Default random number generator object of the type given by $tag$.
..remarks:The random number generator will be default constructed, i.e. with the default seed.
..remarks:This function is NOT thread-safe! Also, data structures using such global state are not thread-safe! Data structures using global random number generator state should use pointers or @Class.Holder|Holder instances@. This way, the random number generator state to be used can be set to be thread-local.
..see:Metafunction.GetDefaultRng
 */

template <typename TRng>
inline TRng &
defaultRng()
{
    static TRng x;
    return x;
}


template <typename T>
inline typename GetDefaultRng<T>::Type &
defaultRng(T const &)
{
    typedef typename GetDefaultRng<T>::Type TRng;
    return defaultRng<TRng>();
}

}  // namespace seqan

#endif  // SEQAN_RANDOM_RANDOM_BASE_H_
