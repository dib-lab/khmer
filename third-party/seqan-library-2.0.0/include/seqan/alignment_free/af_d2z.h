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
// Author: Jonathan Goeke <goeke@molgen.mpg.de>
// ==========================================================================
// This header contains the implementation of the D2z score for
// alignment free sequence comparison.
//
// See Kantorovitz et al. Bioinformatics 2007, Volume23, Issue13,
// Pp. i249-i255.
//
// These functions can be called with alignmentFreeComparison().
// ==========================================================================

// TODO(goeke): const could be added below for the input variables but the function value() in matrix_base (align) is not defined for const. Similarly, the function emittedProbabilty is not defined for const in  statistics_markov_model.h

#ifndef SEQAN_INCLUDE_SEQAN_ALIGNMENT_FREE_AF_D2Z_H_
#define SEQAN_INCLUDE_SEQAN_ALIGNMENT_FREE_AF_D2Z_H_

namespace seqan {

/*
 * _alignmentFreeComparison is called by alignmentFreeComparison() (see alignment_free_comparison.h)
 */
template <typename TStringSet, typename TValue>
void _alignmentFreeComparison(Matrix<TValue, 2> & scoreMatrix,
                              TStringSet const & sequenceSet,
                              AFScore<D2z> const & score)
{
    typedef typename Value<TStringSet>::Type                                    TString;
    typedef typename Value<TString>::Type                                       TAlphabet;
    typedef typename UnmaskedAlphabet_<TAlphabet>::Type                         TUnmaskedAlphabet;
    typedef typename Iterator<TStringSet const>::Type                           TIteratorSet;
    typedef typename Iterator<StringSet<String<unsigned> > >::Type              TIteratorSetUnsigned;
    typedef typename Iterator<StringSet<String<double> > >::Type                TIteratorSetDouble;
    typedef typename Iterator<String<MarkovModel<TUnmaskedAlphabet> > >::Type   TIteratorMarkovModel;

    unsigned seqNumber = length(sequenceSet);

    // Resize scoreMatrix
    setLength(scoreMatrix, 0, seqNumber);
    setLength(scoreMatrix, 1, seqNumber);
    resize(scoreMatrix, (TValue) 0);

    StringSet<String<unsigned> > kmerCounts;
    resize(kmerCounts, seqNumber);

    // Note that there is some code below that looks like copy-and-paste.  However, pulling this out into another
    // function is the only way to get rid of the duplicate lines since we use different types.  After some discussion,
    // weese, goeke and holtgrew agreed that it is probably easier to read and maintain this way than to spread the code
    // over to one more function.
    if (score.bgModelOrder == 0)
    {
        // --------------------------------------------------------------------
        // Order 0 Background Model
        // --------------------------------------------------------------------

        StringSet<String<double> > backgroundFrequencies;
        resize(backgroundFrequencies, seqNumber);

        // Count all kmers and all background nucleotide frequencies and store them in stringSets
        TIteratorSetUnsigned itKmerCounts = begin(kmerCounts);
        TIteratorSetDouble itBackgroundFrequencies = begin(backgroundFrequencies);

        TIteratorSet itSeqSet = begin(sequenceSet);
        for (; itSeqSet < end(sequenceSet); ++itSeqSet)
        {
            countKmers(*itKmerCounts, *itBackgroundFrequencies, *itSeqSet, score.kmerSize);
            ++itKmerCounts;
            ++itBackgroundFrequencies;
        }
        if(score.verbose)
        {
            std::cout << "\ncounted words";
        }
        // Calculate all pairwise scores and store them in scoreMatrix
        for (unsigned rowIndex = 0; rowIndex < seqNumber; ++rowIndex)
        {
            if(score.verbose)
            {
                std::cout << "\nSequence number " << rowIndex;
            }
            for (unsigned colIndex = rowIndex; colIndex < seqNumber; ++colIndex)
            {
                _alignmentFreeCompareCounts(value(scoreMatrix, rowIndex, colIndex), kmerCounts[rowIndex],
                                            backgroundFrequencies[rowIndex], kmerCounts[colIndex],
                                            backgroundFrequencies[colIndex], score);
                value(scoreMatrix, colIndex, rowIndex) = value(scoreMatrix, rowIndex, colIndex);  // Copy symmetric entries
            }
        }
    }
    else
    {
        // --------------------------------------------------------------------
        // Higher Order Background Model
        // --------------------------------------------------------------------

        String<MarkovModel<TUnmaskedAlphabet> > backgroundModels;
        resize(backgroundModels, seqNumber, MarkovModel<TUnmaskedAlphabet>(score.bgModelOrder));
        TIteratorMarkovModel itMM = begin(backgroundModels);
        TIteratorSet itSeqSet = begin(sequenceSet);
        // Count all kmers and all background nucleotide frequencies and store them in StringSets
        TIteratorSetUnsigned itKmerCounts = begin(kmerCounts);

        for (; itSeqSet < end(sequenceSet); ++itSeqSet)
        {
            countKmers(*itKmerCounts, *itMM, *itSeqSet, score.kmerSize);
            ++itKmerCounts;
            if (itMM < end(backgroundModels))
            {
                itMM->_computeAuxiliaryMatrices();
                ++itMM;
            }

        }
        if(score.verbose)
        {
            std::cout << "\ncounted words";
        }
        // Calculate all pairwise scores and store them in scoreMatrix
        for (unsigned rowIndex = 0; rowIndex < seqNumber; ++rowIndex)
        {
            if(score.verbose)
            {
                std::cout << "\nSequence number " << rowIndex;
            }
            for (unsigned colIndex = rowIndex; colIndex < seqNumber; ++colIndex)
            {
                _alignmentFreeCompareCounts(value(scoreMatrix, rowIndex, colIndex), kmerCounts[rowIndex],
                                            backgroundModels[rowIndex], kmerCounts[colIndex],
                                            backgroundModels[colIndex], score);
                value(scoreMatrix, colIndex, rowIndex) = value(scoreMatrix, rowIndex, colIndex);  // Copy symmetric entries
            }
        }
    }
}

/*
 * computeExpectationD2 calculates the expected value of the D2 score given a Bernoulli model,
 * see paper referenced above
 */
template <typename TValue>
double computeExpectationD2(int const len1, int const len2, unsigned const k, TValue const * q1, TValue const * q2)
{
    TValue p2 = 0;
    for (int i = 0; i < 4; i++)
        p2 += q1[i] * q2[i];

    int nbar1 = len1 - k + 1;
    int nbar2 = len2 - k + 1;

    TValue retval = nbar1;
    retval *= nbar2;
    retval *= pow((double)p2, (int)k);
    return retval;
}

/*
 * computeExpectationD2 calculates the expected value of the D2 score given a Markov model,
 * see paper referenced above
 */
template <typename TAlphabet, typename TValue, typename TSpec>
double computeExpectationD2(int const slen1,
                            int const slen2,
                            unsigned const k,
                            MarkovModel<TAlphabet, TValue, TSpec> /*const*/ & bkg1,
                            MarkovModel<TAlphabet, TValue, TSpec> /*const*/ & bkg2,
                            TValue & indicatorexpectation)
{
    unsigned mo = bkg1.order;
    if (mo >= k)
    {
        // Error: Can't suppport markov order greater or equal to word length
        exit(1);
    }

    long emo  = 1 << (2 * mo);  // This is equal to pow(4, mo)
    long ekmo = 1 << (2 * (k - mo));  // This is equal to pow(4, k - mo)

    TValue mean = 0;
    for (long i = 0; i < emo; i++)
    {
        TValue term = value(bkg1.stationaryDistribution, i) * value(bkg2.stationaryDistribution, i);
        TValue subterm = 0;

        for (long j = 0; j < ekmo; j++)
        {
            TValue subprob1 = _computeWordProbGivenPrefix(i, j, bkg1, k, mo);
            TValue subprob2 = _computeWordProbGivenPrefix(i, j, bkg2, k, mo);
            subterm += subprob1 * subprob2;
        }
        mean += term * subterm;
    }

    indicatorexpectation = mean;  // This is E[Y_ij], see paper referenced above.
    int nbar1 = slen1 - k + 1;  // Number of kmers in sequence1
    int nbar2 = slen2 - k + 1;  // Number of kmers in sequence2
    return mean * ((TValue) nbar1 * nbar2);
}

/*
 * computeVarianceD2 calculates the variance of the D2 score given a Bernoulli model,
 * see paper referenced above
 */
template <typename TValue>
double computeVarianceD2(int len1, int len2, unsigned k, TValue * q1, TValue * q2)
{
    int nbar1 = len1 - k + 1;  // Number of kmers in sequence 1
    int nbar2 = len2 - k + 1;  // Number of kmers in sequence 2

    int qbar1 = len1 - 2 * k + 2;  // Number of overlapping kmers in sequence 1
    int qbar2 = len2 - 2 * k + 2;

    TValue p2 = 0, p31 = 0, p32 = 0;
    for (int i = 0; i < 4; i++)
    {
        p2 += q1[i] * q2[i];
        p31 += q1[i] * q2[i] * q1[i];
        p32 += q1[i] * q2[i] * q2[i];
    }

    TValue variance = 0;
    // 'Crabgrass' with l = 0 (= complete overlap), see paper referenced above
    TValue power1 = pow((double)p32, (int)k) - pow((double)p2, 2 * (int)k);
    power1 *= TValue(nbar1) * TValue(qbar2) * TValue(qbar2 - 1);
    TValue power2 = pow((double)p31, (int)k) - pow((double)p2, 2 * (int)k);
    power2 *= TValue(nbar2) * TValue(qbar1) * TValue(qbar1 - 1);
    variance += power1 + power2;

    // 'Crabgrasses' with l > 0, see paper referenced above
    for (unsigned l = 1; l <= k - 1; l++)
        variance += 2 * TValue(nbar1 - l) * TValue(qbar2) * TValue(qbar2 - 1) * (pow((double)p2, (int)(2 * l)) * pow((double)p32, (int)(k - l)) - pow((double)p2, (int)(2 * k))) + 2 * TValue(nbar2 - l) * TValue(qbar1) * TValue(qbar1 - 1) * (pow((double)p2, (int)(2 * l)) * pow((double)p31, (int)(k - l)) - pow((double)p2, (int)(2 * k)));

    // Accordion main diagonal, see paper referenced above
    variance += TValue(nbar1) * TValue(nbar2) * (pow((double)p2, (int)(k)) - pow((double)p2, (int)(2 * k)));
    for (unsigned l = 1; l <= k - 1; l++)
        variance += 2 * TValue(nbar1 - l) * TValue(nbar2 - l) * (pow((double)p2, (int)(k + l)) - pow((double)p2, (int)(2 * k)));

    return variance;
}


/*
 * computeVarianceD2 calculates the variance of the D2 score given a Markov model,
 * see paper referenced above
 */
template <typename TAlphabet, typename TValue, typename TSpec>
double computeVarianceD2(int const slen1,
                         int const slen2,
                         unsigned const k,
                         MarkovModel<TAlphabet, TValue, TSpec> /*const*/ & bkg1,
                         MarkovModel<TAlphabet, TValue, TSpec> /*const*/ & bkg2,
                         TValue indicatorexpectation)
{
    unsigned mo = bkg1.order;
    if (mo >= k)
    {
        // Error: Can't support markov order greater or equal to word length
        exit(1);
    }

    long emo  = 1 << (2 * mo);        // Equivalent to pow(4, mo);
    long ekmo = 1 << (2 * (k - mo));  // Equivalent to pow(4, k - mo);

    int q1 = slen1 - 2 * k + 2;
    int q2 = slen2 - 2 * k + 2;
    int nbar1 = slen1 - k + 1;
    int nbar2 = slen2 - k + 1;

    // Create matrix for Sbctilde
    Matrix<TValue, 2> Sbctilde1;
    setLength(Sbctilde1, 0, emo);
    setLength(Sbctilde1, 1, emo);
    resize(Sbctilde1);

    Matrix<TValue, 2> Sbctilde2;
    setLength(Sbctilde2, 0, emo);
    setLength(Sbctilde2, 1, emo);
    resize(Sbctilde2);

    for (long i = 0; i < emo; i++)
    {
        for (long j = 0; j < emo; j++)
        {
            value(Sbctilde1, i, j) = (q1 * (q1 + 1) / 2) * value(bkg1.stationaryDistribution, j) - (q1 - 1) * value(bkg1._qppp, i, j) - value(bkg1._qppqpp, i, j);
            value(Sbctilde2, i, j) = (q2 * (q2 + 1) / 2) * value(bkg2.stationaryDistribution, j) - (q2 - 1) * value(bkg2._qppp, i, j) - value(bkg2._qppqpp, i, j);
        }
    }

    // Compute sums of word probabilities and star probabilities of x-mers (x between mo + 1 and k) conditional on last or first (resp) mo-word, see paper referenced above
    Matrix<TValue, 2> sump;  // sump[x][wk] is the total word probability of every x-word ending with wk (which is of length mo)
    setLength(sump, 0, k + 1);
    setLength(sump, 1, emo);
    resize(sump, 0.0);

    Matrix<TValue, 2> sumpstar;
    setLength(sumpstar, 0, k + 1);
    setLength(sumpstar, 1, emo);
    resize(sumpstar, 0.0);


    for (unsigned x = 0; x <= k; x++)
    {
        if (x <= mo)
        {

            continue;
        }

        for (long wk = 0; wk < emo; wk++)
        {
            TValue sum = 0;
            long exmo = 1 << (2 * (x - mo));  // Equivalent to pow(4, x - mo);
            for (long wpre = 0; wpre < exmo; wpre++)
            {
                sum += _computeWordProb((wpre << (2 * mo)) + wk, bkg1, x, mo) * _computeWordProb((wpre << (2 * mo)) + wk, bkg2, x, mo);
            }
            value(sump, x, wk) = sum;
        }

        for (long u1 = 0; u1 < emo; u1++)
        {
            TValue sum = 0;
            long exmo = 1 << (2 * (x - mo));  // Equivalent to pow(4, x - mo);
            for (long usuf = 0; usuf < exmo; usuf++)
            {
                sum += _computeWordProbGivenPrefix(u1, usuf, bkg1, x, mo) * _computeWordProbGivenPrefix(u1, usuf, bkg2, x, mo);
            }
            value(sumpstar, x, u1) = sum;
        }
    }

    // Handle the non overlap terms first
    TValue covnonoverlap = 0;
    for (long wk = 0; wk < emo; wk++)
    {
        for (long u1 = 0; u1 < emo; u1++)
        {
            TValue term = value(Sbctilde1, wk, u1) * value(Sbctilde2, wk, u1) * value(sump, k, wk) * value(sumpstar, k, u1);
            covnonoverlap += 4 * term;
        }
    }
    TValue subtractFromNonOverlap = (TValue)q1 * (q1 - 1) * q2 * (q2 - 1);
    subtractFromNonOverlap *= pow(indicatorexpectation, 2);
    covnonoverlap -= subtractFromNonOverlap;
    // Compute 'crabgrass' terms, see paper referenced above
    TValue covcrabgrass = 0;

    // Case 1: overlap >= mo
    for (unsigned m = 1; m <= k - mo; m++)
    {
        long ekm = 1 << (2 * (k - m)); // Equivalent to pow(4, k - m)
        long moonesflag = (1 << (2 * mo)) - 1;
        long kmmoonesflag = (1 << (2 * (k - m - mo))) - 1;
        for (long v = 0; v < ekm; v++)
        {
            long vsuf = v & moonesflag;           // Last morder chars of v
            long vpre = v >> (2 * (k - m - mo));  // First morder chars of v, equivalent to v / pow(4, k - m - mo)
            long vsuf2 = v & kmmoonesflag;        // Remaining k-m-morder chars of v; equivalent to v % pow(4, k - m - mo)
            // In the following, if m = k - mo, _computeWordProbGivenPrefix(.,.,.,k - m, mo) will return 1, as it should
            // Overlap in A, separate in B
            TValue term1 = value(Sbctilde2, vsuf, vpre);
            term1 *= _computeWordProbGivenPrefix(vpre, vsuf2, bkg1, k - m, mo);
            term1 *= pow(_computeWordProbGivenPrefix(vpre, vsuf2, bkg2, k - m, mo), 2);
            term1 *= value(sump, m + mo, vpre);
            term1 *= value(sumpstar, m + mo, vsuf);
            // Overlap in B, separate in A
            TValue term2 = value(Sbctilde1, vsuf, vpre);
            term2 *= _computeWordProbGivenPrefix(vpre, vsuf2, bkg2, k - m, mo);
            term2 *= pow(_computeWordProbGivenPrefix(vpre, vsuf2, bkg1, k - m, mo), 2);
            term2 *= value(sump, m + mo, vpre);
            term2 *= value(sumpstar, m + mo, vsuf);

            covcrabgrass += 4 * q1 * term1 + 4 * q2 * term2;
        }
    }
    // Case 2: overlap < morder
    for (unsigned m = k - mo + 1; m <= k - 1; m++)
    {
        int tlen = 2 * mo - (k - m);
        long et = 1 << (2 * tlen);  // Equivalent to pow(4, tlen);
        long moonesflag = (1 << (2 * mo)) - 1;
        long tmoonesflag = (1 << (2 * (tlen - mo))) - 1;
        for (long t = 0; t < et; t++)
        {
            long tsuf = t & moonesflag;          // Last morder chars of t; equivalent to t % emo;
            long tpre = t >> (2 * (tlen - mo));  // First morder chars of t, equivalent to t / pow(4, tlen - mo)
            long tsuf2 = t & tmoonesflag;        // Remaining tlen-morder chars of t, equivalent to t % etmo;
            // Overlap in A, separate in B
            TValue term1 = value(Sbctilde2, tpre, tsuf);
            term1 *= _computeWordProbGivenPrefix(tpre, tsuf2, bkg1, tlen, mo);
            term1 *= value(sump, k, tpre);
            term1 *= value(sumpstar, k, tsuf);
            // Overlap in B, separate in A
            TValue term2 = value(Sbctilde1, tpre, tsuf);
            term2 *= _computeWordProbGivenPrefix(tpre, tsuf2, bkg2, tlen, mo);
            term2 *= value(sump, k, tpre);
            term2 *= value(sumpstar, k, tsuf);

            covcrabgrass += 4 * q1 * term1 + 4 * q2 * term2;
        }
    }

    // Case 3: m=0 (complete overlap)
    long moonesflag = (1 << (2 * mo)) - 1;
    for (long wpre = 0; wpre < emo; wpre++)
    {
        // First morder chars of w
        TValue term1 = 0;
        TValue term2 = 0;
        for (long wsuf2 = 0; wsuf2 < ekmo; wsuf2++)
        {
            // Remaining k-morder chars of w
            long wsuf = ((wpre << (2 * (k - mo))) + wsuf2) & moonesflag;  // Last morder chars of w
            term1 += 2 * nbar1 * value(Sbctilde2, wsuf, wpre) * _computeWordProbGivenPrefix(wpre, wsuf2, bkg1, k, mo) * pow(_computeWordProbGivenPrefix(wpre, wsuf2, bkg2, k, mo), 2);
            term2 += 2 * nbar2 * value(Sbctilde1, wsuf, wpre) * _computeWordProbGivenPrefix(wpre, wsuf2, bkg2, k, mo) * pow(_computeWordProbGivenPrefix(wpre, wsuf2, bkg1, k, mo), 2);
        }
        covcrabgrass += (term1 + term2) * value(bkg1.stationaryDistribution, wpre) * value(bkg2.stationaryDistribution, wpre);
    }

    // Ignored edge effects, see paper referened above
    // Subtract expectations
    TValue subtractfromcrabgrass = (((2 * q1 * k - nbar1) * (double)q2 * double(q2 - 1)) + ((2 * q2 * k - nbar2) * (double)q1 * double(q1 - 1)));
    subtractfromcrabgrass *= pow(indicatorexpectation, 2);
    covcrabgrass -= subtractfromcrabgrass;

    // Compute the accordion main diagonal terms
    TValue covaccordiondiag = 0;

    for (unsigned m = 1; m <= k - 1; m++)
    {
        long tlen;
        if (m <= mo - 1)
            tlen = m + mo;
        else
            tlen = 2 * mo - 1;
        long et = 1 << (2 * tlen);  // Equivalent to pow(4, tlen)
        long moonesflag = (1 << (2 * mo)) - 1;
        long tmoonesflag = (1 << (2 * (tlen - mo))) - 1;

        TValue mterm = 0;
        for (long t = 0; t < et; t++)
        {
            TValue term = 1;
            long tpre = t >> (2 * (tlen - mo));  // First morder chars of t, equivalent to t / pow(4, tlen - mo)
            long tsuf2 = t & tmoonesflag;        // The remaining tlen - morder chars of t, equivalent to t % etmo;
            term *= _computeWordProbGivenPrefix(tpre, tsuf2, bkg1, tlen, mo) * _computeWordProbGivenPrefix(tpre, tsuf2, bkg2, tlen, mo);
            term *= value(sump, k, tpre);
            if (m >= mo)
            {
                long tsuf = t & moonesflag;  // Last morder chars of t, equivalent to t%emo;
                term *= value(sumpstar, m + 1, tsuf);
            }
            mterm += term;
        }
        covaccordiondiag += 2 * (nbar1 - m) * (nbar2 - m) * mterm;
    }
    covaccordiondiag += nbar1 * nbar2 * indicatorexpectation;  // Complete overlap term
    covaccordiondiag -= (nbar1 * nbar2 * (2 * k - 1) - (nbar1 + nbar2) * k * (k - 1) + (k - 1) * k * (2 * k - 1) / 3) * pow(indicatorexpectation, 2);
    TValue variance = covnonoverlap + covcrabgrass + covaccordiondiag;
    return variance;

}

/*
 * Calculate pairwise score given the counts of all kmers and the background Bernoulli models
 */
template <typename TValue, typename TStringBG>
void _alignmentFreeCompareCounts(TValue & result,
                                String<unsigned> const & kmerCounts1,
                                TStringBG const & backgroundFrequencies1,
                                String<unsigned> const & kmerCounts2,
                                TStringBG const & backgroundFrequencies2,
                                AFScore<D2z> const & score)
{
    typedef typename Value<TStringBG>::Type TValueBG;
    TValueBG sum = 0;
    unsigned len1 = score.kmerSize - 1;
    unsigned len2 = score.kmerSize - 1;
    unsigned nvals = length(kmerCounts1);

    for (unsigned l = 0; l < nvals; l++)
    {
        len1 += kmerCounts1[l];
        len2 += kmerCounts2[l];

        sum += kmerCounts1[l] * kmerCounts2[l];
    }

    TValueBG q1[4];
    TValueBG q2[4];
    for (int l = 0; l < 4; l++)
    {
        q1[l] = backgroundFrequencies1[l];
        q2[l] = backgroundFrequencies2[l];
    }

    // Compute expected value and variance (IID)
    double E = computeExpectationD2(len1, len2, score.kmerSize, q1, q2);
    double var = computeVarianceD2(len1, len2, score.kmerSize, q1, q2);

    if ((var <= 0))
    {
        if(score.verbose)
        {
            std::cout << "Error: negative variance\n";
        }
        result = 0;
        return;
    }
    // Calculate z-score
    result = (TValue) (sum - E) / pow(var, 0.5);
}

/*
 * Calculate pairwise score given the counts of all kmers and the background Markov models
 */
template <typename TAlphabet, typename TValue, typename TSpec>
void _alignmentFreeCompareCounts(TValue & result,
                                String<unsigned> const & kmerCounts1,
                                MarkovModel<TAlphabet, TValue, TSpec> /*const*/ & bgModel1,
                                String<unsigned> const & kmerCounts2,
                                MarkovModel<TAlphabet, TValue, TSpec> /*const*/ & bgModel2,
                                AFScore<D2z> const & score)
{
    unsigned nvals = length(kmerCounts1);
    int sum = 0;
    int sumCounts1 = score.kmerSize - 1;
    int sumCounts2 = score.kmerSize - 1;

    for (unsigned l = 0; l < nvals; l++)
    {
        sumCounts1 += kmerCounts1[l];
        sumCounts2 += kmerCounts2[l];
        sum += value(kmerCounts1, l) * value(kmerCounts2, l);  // Calculate the inner product
    }

    TValue D2 = (TValue) sum;

    // Compute mean and variance
    TValue indicatorexpectation = 0;

    double E = computeExpectationD2(sumCounts1, sumCounts2, score.kmerSize, bgModel1, bgModel2, indicatorexpectation);
    double var = computeVarianceD2(sumCounts1, sumCounts2, score.kmerSize, bgModel1, bgModel2, indicatorexpectation);

    if (var <= 0)
    {
        result = 0;
        return;
    }
    // Return z-score of D2
    result = (D2 - E) / pow(var, 0.5);
}

/*
 * Compute the word probability given a Markov model
 * see paper referenced on top
 */
template <typename TAlphabet, typename TValue, typename TSpec>
double _computeWordProb(long const word, MarkovModel<TAlphabet, TValue, TSpec> /*const*/ & bkg, unsigned const k, int const mo)
{
    // mo needs to be passed, since that decides the prefix
    // Similarly, this function needs to know the total length of the word k
    long prefix = word >> (2 * (k - mo));  // Equivalent to word / pow(4, k - mo); gets the first mo chars of word
    long suffix = word & ((1 << (2 * (k - mo))) - 1);  // Get the last (k - mo) chars of word
    return value(bkg.stationaryDistribution, prefix) * _computeWordProbGivenPrefix(prefix, suffix, bkg, k, mo);
}

/*
 * Compute the word probability given a Markov model
 * see paper referenced on top
 */
template <typename TAlphabet, typename TValue, typename TSpec>
double _computeWordProbGivenPrefix(long const prefix, long const suffix, MarkovModel<TAlphabet, TValue, TSpec> /*const*/ & bkg, unsigned const k, unsigned const mo)
{
    // Computes p_ * (prefix suffix)
    // Calculate only if k - mo >= 1; otherwise return 1
    if (k - mo <= 0)
        return 1;

    TValue prob = 1;
    long pre = prefix;
    long jcopy = suffix;

    // Iterate through successive positions of the suffix
    for (unsigned l = 0; l < k - mo; l++)
    {
        // Get the first char of jcopy, which is of length k - mo; equivalent to jcopy / pow(4, k - mo - 1)
        long jcopyfirst = jcopy >> (2 * (k - mo - 1));
        long suf = ((pre & ((1 << (2 * (mo - 1))) - 1)) << 2) + jcopyfirst;

        // prob *=P[pre][suf]; P is the transition matrix
        prob *= value(bkg.transition, pre, suf);
        // Erase the first char of jcopy, and preserve its length at k - mo by shifting one char to the left; equivalent to (jcopy % pow(4, k - mo - 1)) * 4;
        jcopy = (jcopy & ((1 << (2 * (k - mo - 1))) - 1)) << 2;
        pre = suf;
    }
    return prob;
}

}  // namespace seqan

#endif  // SEQAN_INCLUDE_SEQAN_ALIGNMENT_FREE_AF_D2Z_H_
