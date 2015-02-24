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
// Code for geometrically distributed random number generation, p=0.5.
// ==========================================================================

#ifndef SEQAN_RANDOM_RANDOM_GEOMETRIC_H_
#define SEQAN_RANDOM_RANDOM_GEOMETRIC_H_

namespace seqan {

// ===========================================================================
// Forwards, Tags.
// ===========================================================================

// Specialization Tag for geometric distribution.
struct GeometricFairCoin {};

// ===========================================================================
// Classes
// ===========================================================================

/*!
 * @class GeometricPdf
 * @headerfile <seqan/random.h>
 * @extends Pdf
 * @signature template <>
 *            class Pdf<Geometric>;
 * @brief Geometric probability density function with <i>p = 0.5</i>.
 *
 * The PDF is implemented efficiently without using any floating point arithmetics, just some bit operations.
 *
 * @fn GeometricPdf::Pdf
 * @brief Constructor
 *
 * @signature Pdf::Pdf();
 */

template <>
class Pdf<GeometricFairCoin>
{
public:

    Pdf() { SEQAN_CHECKPOINT; }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

template <>
struct Value<Pdf<GeometricFairCoin> >
{
    typedef unsigned Type;
};

template <>
struct Value<const Pdf<GeometricFairCoin> > : Value<Pdf<GeometricFairCoin> > {};

// ===========================================================================
// Functions
// ===========================================================================

template <typename TRNG>
inline
typename Value<Pdf<GeometricFairCoin> >::Type
pickRandomNumber(TRNG & rng, Pdf<GeometricFairCoin> const & /*pdf*/)
{
    SEQAN_CHECKPOINT;

    const int RG_IB1 = 1;
    const int RG_IB2 = 2;
    const int RG_IB5 = 16;
    const int RG_IB18 = 131072;
    const int RG_MASK = RG_IB1 + RG_IB2 + RG_IB5;

    typename Value<TRNG>::Type seed = pickRandomNumber(rng);
    typename Value<Pdf<GeometricFairCoin> >::Type value = 0;

    while (true) {
        if ((seed & RG_IB18)) {
            seed = ((seed ^ RG_MASK) << 1) | RG_IB1;
            ++value;
        } else {
            seed <<= 1;
            break;
        }
    }

    return value;
}

}  // namespace seqan

#endif  // SEQAN_RANDOM_RANDOM_GEOMETRIC_H_
