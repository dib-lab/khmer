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
// Code for random number generation following beta distribution.
// ==========================================================================

#ifndef SEQAN_RANDOM_RANDOM_BETA_H_
#define SEQAN_RANDOM_RANDOM_BETA_H_

#include "random_beta_kfunc.h"

namespace seqan {

// ===========================================================================
// Forwards, Tags.
// ===========================================================================

// Specialization Tag for normal distribution.
struct Beta_;
typedef Tag<Beta_> Beta;

// Selection of construction with alpha and beta.
struct AlphaBeta_;
typedef Tag<AlphaBeta_> AlphaBeta;

// Selection of construction with mu and sigma.
struct MeanStdDev_;
typedef Tag<MeanStdDev_> MeanStdDev;

// ===========================================================================
// Classes
// ===========================================================================

/*!
 * @class BetaPdf Beta Pdf
 * @extends Pdf
 * @headerfile <seqan/random.h>
 * @brief Probability density function for the Beta function.
 *
 * @signature template <>
 *            class Pdf<Beta>;
 *
 *
 * @fn BetaPdf::Pdf
 * @brief Constructor for the beta Pdf.
 *
 * Use the tags <tt>AlphaBeta</tt> and <tt>MeanStdDev</tt> to select the meaning of the two parameters.
 *
 * @signature Pdf::Pdf(alpha, beta[,  AlphaBeta()]);
 * @signature Pdf::Pdf(mu,    sigma, MeanStdDev());
 *
 * @param[in] alpha Alpha value for the beta distribution, <tt>double</tt>.
 * @param[in] beta  Beta value for the beta distribution, <tt>double</tt>.
 * @param[in] mu    Mean value for the beta distribution, <tt>double</tt>.
 * @param[in] sigma Standard deviation value for the beta distribution, <tt>double</tt>.
 */

template <>
class Pdf<Beta>
{
public:
    double _alpha;
    double _beta;

    Pdf(double mu, double sigma, MeanStdDev const & /*tag*/)
            : _alpha(((1 - mu) / sigma / sigma - 1 / mu) * mu * mu),
              _beta(_alpha * (1 / mu - 1))
    {}

    Pdf(double alpha, double beta, AlphaBeta const & /*tag*/) :
            _alpha(alpha), _beta(beta)
    {}

    Pdf(double alpha, double beta) :
            _alpha(alpha), _beta(beta)
    {}
};

// ===========================================================================
// Metafunctions
// ===========================================================================

template <>
struct Value<Pdf<Beta> >
{
    typedef double Type;
};

template <>
struct Value<const Pdf<Beta> > : Value<Pdf<Beta> > {};

// ===========================================================================
// Functions
// ===========================================================================

template <typename TRNG>
inline
typename Value<Pdf<Beta> >::Type
pickRandomNumber(TRNG & rng, Pdf<Beta> const & pdf)
{
    Pdf<Uniform<double> > pdfUniform(0, 1);
    return 1 - kf_betai(pdf._alpha, pdf._beta, pickRandomNumber(rng, pdfUniform));
}

}  // namespace seqan

#endif  // SEQAN_RANDOM_RANDOM_BETA_H_
