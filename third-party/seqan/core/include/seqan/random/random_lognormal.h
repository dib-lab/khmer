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
// Code for log-normally distributed random number generation.
// ==========================================================================

#ifndef SEQAN_RANDOM_RANDOM_LOGNORMAL_H_
#define SEQAN_RANDOM_RANDOM_LOGNORMAL_H_

namespace seqan {

// ===========================================================================
// Forwards, Tags.
// ===========================================================================

// Forward-declarations.
struct Normal_;
typedef Tag<Normal_> Normal;

// Specialization Tag for log-normal distribution.
struct LogNormal_;
typedef Tag<LogNormal_> LogNormal;

/**
.Tag.Lognormal Construction:
..summary:Specify parameters for lognormal distribution construction.
..cat:Random
..include:seqan/random.h
..tag.MuSigma:
...summary:Tag to specify that the given parameters are mu and sigma of the underlying normal distribution for lognormal distributions.
..tag.MeanStdDev:
...summary:Tag to specify that the given parameters are mean an standard deviation of the lognormal distribution.
..wiki:Tutorial/Randomness|Tutorial: Randomness
..see:Spec.Log-Normal Pdf
*/

struct MuSigma_;
typedef Tag<MuSigma_> MuSigma;
struct MeanStdDev_;
typedef Tag<MeanStdDev_> MeanStdDev;

// ===========================================================================
// Classes
// ===========================================================================

/**
.Spec.Log-Normal Pdf
..general:Class.Pdf
..summary:Log-normal probability density function.
..remark:Note that you can construct this either with mu/sigma of the underlying normal distribution (default) or with the mean and standard deviation of the log-normal distribution.
..cat:Random
..include:seqan/random.h
..wiki:Tutorial/Randomness|Tutorial: Randomness
*/

template <>
class Pdf<LogNormal>
{
public:
    Pdf<Normal> _normalDist;

/**
.Memfunc.Log-Normal Pdf#Pdf
..class:Spec.Log-Normal Pdf
..summary:Constructor for log-normal Pdf.
Log-normal PDFs can either be initialized by the mean and standard deviation of the underlying normal distribution or directly of the log-normal distribution.
..signature:Pdf<LogNormal>(mu, sigma[, MuSigma()])
..signature:Pdf<LogNormal>(mean, stdDev, MeanStdDev())
..param.mu:Mean of the underlying normal distribution.
...type:nolink:double
..param.sigma:Standard deviation of the underlying normal distribution.
...type:nolink:double
..param.mean:Mean of the log-normal distribution.
...type:nolink:double
..param.stdDev:Standard deviation of the log-normal distribution.
...type:nolink:double
..see:Tag.Lognormal Construction.tag.MuSigma
..see:Tag.Lognormal Construction.tag.MeanStdDev
*/
    Pdf(double mu, double sigma, MuSigma const &)
            : _normalDist(mu, sigma)
    {
        SEQAN_CHECKPOINT;
    }

    Pdf(double mean, double stddev, MeanStdDev const &)
            : _normalDist(::std::log(mean) - 0.5 * ::std::log(1.0 + stddev * stddev / mean / mean),
                          ::std::sqrt(::std::log(1.0 + stddev * stddev / mean / mean)))
    {
        SEQAN_CHECKPOINT;
    }

    Pdf(double mu, double sigma)
            : _normalDist(mu, sigma)
    {
        SEQAN_CHECKPOINT;
    }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

template <>
struct Value<Pdf<LogNormal> >
{
    typedef double Type;
};

template <>
struct Value<const Pdf<LogNormal> > : Value<Pdf<LogNormal> > {};

// ===========================================================================
// Functions
// ===========================================================================

template <typename TRandomNumberGenerator>
inline
typename Value<Pdf<LogNormal> >::Type
pickRandomNumber(TRandomNumberGenerator & rng, Pdf<LogNormal> const & pdf)
{
    SEQAN_CHECKPOINT;
    return exp(pickRandomNumber(rng, pdf._normalDist));
}

}  // namespace seqan

#endif  // SEQAN_RANDOM_RANDOM_LOGNORMAL_H_
