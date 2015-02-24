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
// Functor for random number generation.
// ==========================================================================

#ifndef SEQAN_RANDOM_RANDOM_RNG_FUNCTOR_H_
#define SEQAN_RANDOM_RANDOM_RNG_FUNCTOR_H_

namespace seqan {

// ===========================================================================
// Forwards, Tags.
// ===========================================================================

// Tag for selecting the Rng functor specialization.
template <typename TRng, typename TPdf>
struct RngFunctor {};

// ===========================================================================
// Classes
// ===========================================================================

/*!
 * @class RngFunctor
 * @extends Rng
 * @headerfile <seqan/random.h>
 * @brief Functor wrapper for random number generation.
 *
 * @signature template <typename TRng, typename TPdf>
 *            class Rng<RngFunctor<TRng, TPdf> >;
 *
 * @tparam TRng The random number generator type to use.
 * @tparam TPdf The probability density function type to use.
 */

template <typename TRng, typename TPdf>
class Rng<RngFunctor<TRng, TPdf> >
{
public:
    TRng & _rng;
    TPdf & _pdf;

/*!
 * @fn RngFunctor::Rng
 * @headerfile <seqan/random.h>
 * @brief Constructor.
 *
 * @signature Rng::Rng(rng, pdf);
 *
 * @param[in] rng A reference to the underlying Rng to use.
 * @param[in] pdf A reference to the underlying Pdf to use.
 */

    Rng(TRng & rng, TPdf & pdf)
        : _rng(rng), _pdf(pdf)
    {}

    inline
    typename Value<TPdf>::Type
    operator()()
    {
        SEQAN_CHECKPOINT;
        return pickRandomNumber(_rng, _pdf);
    }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

template <typename TRng, typename TPdf>
struct Value<RngFunctor<TRng, TPdf> > : Value<TPdf> {};

template <typename TRng, typename TPdf>
struct Value<RngFunctor<TRng, TPdf> const> : Value<TPdf> {};

template <typename TRng, typename TPdf>
struct MaxValue<RngFunctor<TRng, TPdf> > : MaxValue<TPdf> {};

template <typename TRng, typename TPdf>
struct MaxValue<RngFunctor<TRng, TPdf> const> : MaxValue<TPdf> {};

template <typename TRng, typename TPdf>
struct MinValue<RngFunctor<TRng, TPdf> > : MinValue<TPdf> {};

template <typename TRng, typename TPdf>
struct MinValue<RngFunctor<TRng, TPdf> const> : MinValue<TPdf> {};

// ===========================================================================
// Functions
// ===========================================================================

template <typename TRng, typename TPdf>
inline unsigned
pickRandomNumber(Rng<RngFunctor<TRng, TPdf> > & rng)
{
    SEQAN_CHECKPOINT;

    return pickRandomNumber(rng._rng, rng._pdf);
}

}  // namespace seqan

#endif  // SEQAN_RANDOM_RANDOM_RNG_FUNCTOR_H_
