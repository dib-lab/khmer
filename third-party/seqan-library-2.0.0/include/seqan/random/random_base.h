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

/*!
 * @class Rng
 * @headerfile <seqan/random.h>
 * @brief Random Number Generator
 *
 * @signature template <[typename TSpec]>
 *                     class Rng;
 *
 * @tparam TSpec Tag for selecting the specialization.  Defaults to @link MersenneTwisterRng @endlink.
 *
 * @section Examples
 *
 * The following code shows how to generate random numbers and shuffle a text.
 *
 * @include demos/random/random.cpp
 *
 * @code{.console}
 * pickRandomNumber(rng) == 1608637542
 * pickRandomNumber(rng, uniformDouble) == 0.950714
 * pickRandomNumber(rng, uniformInt) == 27
 * pickRandomNumber(rng, normal) == 0.419823
 * pickRandomNumber(rng, logNormal) == 1.22431
 * pickRandomNumber(rng, logNormal2) == 2.78004
 * pickRandomNumber(rng, logNormal3) == 0.00155248
 * shuffle("Hello World!") ==  o!reWlloHld
 * @endcode
 *
 * @see http://seqan.readthedocs.org/en/develop/Tutorial/Randomness.html
 */

template <typename TSpec = MersenneTwister>
class Rng;


// ===========================================================================
// Classes
// ===========================================================================

/*!
 * @fn Rng::operator()
 * @brief Function all operator.
 *
 * @signature TValue Rng::operator();
 *
 * @return TValue Random number, TValue can be retrieved with Rng#Value.
 */

/*!
 * @class Pdf
 * @headerfile <seqan/random.h>
 * @brief Probability Density Function
 *
 * @signature template <typename TSpec>
 *            class Pdf;
 *
 * @tparam TSpec Tag for selecting the specialization.
 *
 * This class is used together with @link Rng @endlink in the function @link Rng#pickRandomNumber @endlink. See the
 * <a href="http://seqan.readthedocs.org/en/develop/Tutorial/Randomness.html">SeqAn Randomness Tutorial</a> for more
 * details.
 */

template <typename TSpec>
class Pdf;


// ===========================================================================
// Metafunctions
// ===========================================================================

/*!
 * @mfn Pdf#Value
 * @brief Value type of a Pdf.
 *
 * @signature Value<TPdf>::Type
 *
 * @tparam TPdf The Pdf for the value type.
 *
 * @return Type The value type of the Pdf.
 */

// specification only

/*!
 * @mfn Rng#Value
 * @brief Value type of a Rng.
 *
 * @signature Value<TRng>::Type
 *
 * @tparam TRng the Rng to get the value type for.
 *
 * @return Type the value type for the Rng.
 */

/*!
 * @mfn Rng#MinValue
 * @brief Smallest value that a Rng can return.
 *
 * @signature MinValue<TRng>::VALUE;
 *
 * @tparam TRng The Rng object to get the smallest value for.
 *
 * @return VALUE The smallest value a Rng can return.
 */

/*!
 * @mfn Rng#MaxValue
 * @brief Largest value that a Rng can return.
 *
 * @signature MaxValue<TRng>::VALUE;
 *
 * @tparam TRng The Rng object to get the largest value for.
 *
 * @return VALUE The largest value a Rng can return.
 */

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

/*!
 * @mfn GetDefaultRng
 * @headerfile <seqan/random.h>
 * @brief Return the default Rng to use in a given class, specialization, or algorithm.
 *
 * @signature GetDefaultRng<T>::Type;
 *
 * @tparam T The type or tag to get the default Rng for.
 *
 * @return Type The Rng type to use.
 *
 * @see defaultRng
 */

template <typename T>
struct GetDefaultRng
{
    typedef Rng<MersenneTwister> Type;
};

// ===========================================================================
// Functions
// ===========================================================================

/*!
 * @fn Rng#pickRandomNumber
 * @brief Pick a random number using a random number generator, possibly using a probability density function.
 *
 * @signature TValue pickRandomNumber(rng[, pdf]);
 *
 * @param[in,out] rng The Rng to use.
 * @param[in]     pdf The probability density function to use.
 *
 * @return TValue A random number.  TValue is the value type of the rng if pdf is not given.  If pdf is given
 *                then it is of the value type of pdf.
 *
 * @section Remarks
 *
 * For more details see the <a href="http://seqan.readthedocs.org/en/develop/Tutorial/Randomness.html">SeqAn
 * Tutorial on Randomness</a>.
 */

// specification only

/*!
 * @fn defaultRng
 * @headerfile <seqan/random.h>
 * @brief Return the default random number generator object of a given type.
 *
 * @signature TRng defaultRng<TRng>();
 *
 * @tparam TRng The random number generator type to construct.
 *
 * @return TRng Default random number generator of the given type TRng.
 *
 * @section Remarks
 *
 * The random number generator will be default constructed, i.e. with the default seed.
 *
 * This function is NOT thread-safe!  Also, data structures and functions using defaultRng are not thread-safe.  Data
 * structures using global random number generator state should use pointers.  This way, the random number generator
 * state to be used can be set to be thread-local.
 *
 * @see GetDefaultRng
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
