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
// Code for normally distributed random number generation.
// ==========================================================================

#ifndef SEQAN_RANDOM_RANDOM_NORMAL_H_
#define SEQAN_RANDOM_RANDOM_NORMAL_H_

namespace seqan {

// ===========================================================================
// Forwards, Tags.
// ===========================================================================

// Specialization Tag for normal distribution.
struct Normal_;
typedef Tag<Normal_> Normal;

// ===========================================================================
// Classes
// ===========================================================================

/**
.Spec.Normal Pdf
..signature:Pdf<Normal>
..general:Class.Pdf
..summary:Normal probability density function.
..cat:Random
..include:seqan/random.h
..wiki:Tutorial/Randomness|Tutorial: Randomness
*/

template <>
class Pdf<Normal>
{
public:
    double _mu;
    double _sigma;

/**
.Memfunc.Normal Pdf#Pdf
..class:Spec.Normal Pdf
..summary:Constructor for normal Pdf.
..signature:Pdf<Normal>(mu, sigma)
..param.mu:Mean of the normal distribution.
...type:nolink:double
..param.sigma:Standard deviation of the normal distribution.
...type:nolink:double
*/
    Pdf(double mu, double sigma)
            : _mu(mu), _sigma(sigma)
    {
        SEQAN_CHECKPOINT;
    }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

template <>
struct Value<Pdf<Normal> >
{
    typedef double Type;
};

template <>
struct Value<const Pdf<Normal> > : Value<Pdf<Normal> > {};

// ===========================================================================
// Functions
// ===========================================================================

static const double SEQAN_NV_MAGICCONST = 1.7155277699214135;  // == 4 * exp(-0.5)/sqrt(2.0)

/*
..summary:Pick a normally distributed random number.
*/
template <typename TRNG>
inline
typename Value<Pdf<Normal> >::Type
pickRandomNumber(TRNG & rng, Pdf<Normal> const & pdf)
{
    SEQAN_CHECKPOINT;

    // Normal Distribution Heuristics, ported from Python.
    //
    // Kinderman and Monahan method. Reference: Kinderman, A.J. and
    // Monahan, J.F., "Computer generation of random variables using
    // the ratio of uniform deviates", ACM Trans Math Software, 3,
    // (1977), pp257-260.

    double z;
    Pdf<Uniform<double> > pdfUniform(0, 1);
    while (true) {
        double u1 = pickRandomNumber(rng, pdfUniform);
        double u2 = 1 - pickRandomNumber(rng, pdfUniform);
        z = SEQAN_NV_MAGICCONST * (u1 - 0.5) / u2;
        double zz = z * z / 4.0;
        if (zz < -::std::log10(u2))
            break;
    }
    return pdf._mu + z * pdf._sigma;
}

}  // namespace seqan

#endif  // SEQAN_RANDOM_RANDOM_NORMAL_H_
