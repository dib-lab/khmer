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
// Code for uniformly distributed random number generation.
// ==========================================================================

#ifndef SEQAN_RANDOM_RANDOM_UNIFORM_H_
#define SEQAN_RANDOM_RANDOM_UNIFORM_H_

namespace seqan {

// ===========================================================================
// Forwards, Tags.
// ===========================================================================

// Specialization tag for uniform distribution.
template <typename T>
struct Uniform;

// ===========================================================================
// Classes
// ===========================================================================

/*!
 * @class UniformPdf Uniform Pdf
 * @headerfile <seqan/random.h>
 * @extends Pdf
 * @brief Uniform distribution probability density function over a closed interval [min, max];
 *
 * @signature template <typename T>
 *            class Pdf<Uniform<T> >;
 *
 * @tparam T The number type to use (integer or floating point).
 */

template <typename T>
class Pdf<Uniform<T> >
{
public:
    T _min;
    T _max;

/*!
 * @fn UniformPdf::Pdf
 * @brief Constructor for uniform Pdf.
 *
 * @signature Pdf::Pdf(min, max);
 *
 * @param[in] min The smallest value of the interval, of the Value type of the Pdf.
 * @param[in] max The largest value of the interval, of the Value type of the Pdf.
 */

    Pdf(T min, T max)
            : _min(min), _max(max)
    {
        SEQAN_CHECKPOINT;
        SEQAN_ASSERT_LEQ(_min, _max);
    }
};

// Specialization for bools do not need min/max.
template <>
class Pdf<Uniform<bool> >
{
public:
    Pdf() {}
};

// ===========================================================================
// Metafunctions
// ===========================================================================

template <typename T>
struct Value<Pdf<Uniform<T> > >
{
    typedef T Type;
};

template <typename T>
struct Value<const Pdf<Uniform<T> > > : Value<Pdf<Uniform<T> > > {};

// ===========================================================================
// Functions
// ===========================================================================

// Pick an integral random number uniformly distributed.
template <typename TRNG, typename T>
inline
typename Value<Pdf<Uniform<T> > >::Type
_pickRandomNumber(TRNG & rng, Pdf<Uniform<T> > const & pdf, True const &)
{
    SEQAN_CHECKPOINT;
    typename Value<TRNG>::Type limit = (MaxValue<TRNG>::VALUE / (pdf._max - pdf._min)) * (pdf._max - pdf._min);
    typename Value<TRNG>::Type x;
    do {
        x = pickRandomNumber(rng);
    } while (x > limit);
    T y = x % (pdf._max - pdf._min + 1);
    return y + pdf._min;
}

// Pick a continuous random number uniformly distributed.
template <typename TRNG, typename T>
inline
typename Value<Pdf<Uniform<T> > >::Type
_pickRandomNumber(TRNG & rng, Pdf<Uniform<T> > const & pdf, False const &)
{
    SEQAN_CHECKPOINT;
    T x = static_cast<T>(pickRandomNumber(rng) - MinValue<TRNG>::VALUE);
    x /= static_cast<T>(MaxValue<TRNG>::VALUE) - static_cast<T>(MinValue<TRNG>::VALUE);
    return pdf._min + x * (pdf._max - pdf._min);
}

template <typename TRNG, typename T>
inline
typename Value<Pdf<Uniform<T> > >::Type
pickRandomNumber(TRNG & rng, Pdf<Uniform<T> > const & pdf)
{
    SEQAN_CHECKPOINT;
    if (pdf._min == pdf._max)
        return pdf._min;
    return _pickRandomNumber(rng, pdf, typename IsInteger<T>::Type());
}

// Specialization for picking a random bool.
template <typename TRNG>
inline
typename Value<Pdf<Uniform<bool> > >::Type
pickRandomNumber(TRNG & rng, Pdf<Uniform<bool> > const &)
{
    SEQAN_CHECKPOINT;
    typename Value<TRNG>::Type x = pickRandomNumber(rng);
    return x % 2;
}

}  // namespace seqan

#endif  // SEQAN_RANDOM_RANDOM_UNIFORM_H_
