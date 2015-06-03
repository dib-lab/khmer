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
// This header contains the implementation of the D2star score for alignment
// free sequence comparison.
//
// See: Reinert et al. J Comput Biol. 2009 Dec;16(12):1615-34.
//
// These functions can be called with alignmentFreeComparison().
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGNMENT_FREE_AF_D2STAR_ORIGINAL_H_
#define SEQAN_INCLUDE_SEQAN_ALIGNMENT_FREE_AF_D2STAR_ORIGINAL_H_

namespace seqan {

/*
 * _alignmentFreeComparison is called by alignmentFreeComparison() (see alignment_free_comparison.h)
 */
template <typename TValue, typename TStringSet>
void _alignmentFreeComparison(Matrix<TValue, 2> & scoreMatrix,
                              TStringSet const & sequenceSet,
                              AFScore<D2Star> const & score)
{

    //typedef typename Iterator<TStringSet const>::Type             TIteratorSet;
    //typedef typename Iterator<StringSet<String<double> > >::Type  TIteratorSetDouble;

    //typedef Matrix<TValue, 2> TMatrix;

    unsigned seqNumber = length(sequenceSet);

    // Resize the scoreMatrix
    setLength(scoreMatrix, 0, seqNumber);
    setLength(scoreMatrix, 1, seqNumber);
    resize(scoreMatrix, (TValue) 0);
    StringSet<String<double> > standardisedKmerCounts;
    resize(standardisedKmerCounts, seqNumber);

    // Calculate all pairwise scores and store them in scoreMatrix
    for (unsigned rowIndex = 0; rowIndex < seqNumber; ++rowIndex)
    {
        if(score.verbose)
        {
            std::cout << "\nSequence number " << rowIndex;
        }
        for (unsigned colIndex = rowIndex; colIndex < seqNumber; ++colIndex)
        {
            _d2star(value(scoreMatrix, rowIndex, colIndex), sequenceSet[rowIndex], sequenceSet[colIndex], score);
            value(scoreMatrix, colIndex, rowIndex) = value(scoreMatrix, rowIndex, colIndex);  // Copy symmetric entries
        }
    }
}

/*
 * _d2star calculates the pairwise score of two sequences according to the paper referenced above.
 */
template <typename TValue, typename TSequence>
void _d2star(TValue & result,
             TSequence const & sequence1,
             TSequence const & sequence2,
             AFScore<D2Star> const & score)
{
    typedef typename Value<TSequence>::Type              TAlphabet;
    typedef typename UnmaskedAlphabet_<TAlphabet>::Type  TUnmaskedAlphabet;

    TValue missing = -pow(10.0, 10);
    TSequence seq1seq2;
    append(seq1seq2, sequence1);
    append(seq1seq2, sequence2);
    result = 0.0;

    // Note that there is some code below that looks like copy-and-paste.  However, pulling this out into another
    // function is the only way to get rid of the duplicate lines since we use different types.  After some discussion,
    // weese, goeke and holtgrew agreed that it is probably easier to read and maintain this way than to spread the code
    // over to one more function.
    if (score.bgModelOrder == 0)
    {
        // --------------------------------------------------------------------
        // Order 0 Background Model
        // --------------------------------------------------------------------

        String<unsigned> kmerCounts1;
        String<unsigned> kmerCounts2;
        String<unsigned> backgroundCounts;
        String<double> backgroundFrequencies;
        resize(backgroundFrequencies, 4, 0);
        countKmers(kmerCounts1, sequence1, score.kmerSize);
        countKmers(kmerCounts2, sequence2, score.kmerSize);
        countKmers(backgroundCounts, seq1seq2, 1);
        int sumBG = 0;
        for (unsigned i = 0; i < length(backgroundCounts); ++i)
        {
            sumBG += backgroundCounts[i];
        }
        for (unsigned i = 0; i < length(backgroundCounts); ++i)
        {
            backgroundFrequencies[i] = backgroundCounts[i] / ((double)sumBG);
        }
        unsigned nvals = length(kmerCounts1);  // Number of kmers
        int len1 = 0;
        int len2 = 0;

        for (unsigned l = 0; l < nvals; l++)
        {
            len1 += kmerCounts1[l];
            len2 += kmerCounts2[l];
        }

        String<TValue> probabilities;  // String of TValue to store the word probabilities p_w
        resize(probabilities, nvals, missing);

        for (unsigned i = 0; i < nvals; ++i)
        {
            TValue p_w = 1;  // Probability of kmer

            String<TUnmaskedAlphabet> w;
            unhash(w, i, score.kmerSize);
            calculateProbability(p_w, w, backgroundFrequencies);
            TValue variance1 = 0.0;
            TValue variance2 = 0.0;

            variance1 = pow(len1 * p_w, 0.5);
            variance2 = pow(len2 * p_w, 0.5);

            // Test if variance is larer than 0 and smaller than inf before dividing
            if ((variance1 > missing) && (variance1 < pow(10.0, 10)))
            {
                if (p_w > 0)
                {
                    TValue stCount1 = (kmerCounts1[i] - p_w * len1) / variance1;
                    TValue stCount2 = (kmerCounts2[i] - p_w * len2) / variance2;
                    result += stCount1 * stCount2;
                }
            }
        }
    }
    else
    {
        // --------------------------------------------------------------------
        // Higher Order Background Model
        // --------------------------------------------------------------------

        String<unsigned> kmerCounts1;
        String<unsigned> kmerCounts2;
        StringSet<String<TUnmaskedAlphabet> > bgSequences;
        stringToStringSet(bgSequences, seq1seq2);  // Create unmasked sequences
        MarkovModel<TUnmaskedAlphabet, TValue> backgroundModel(score.bgModelOrder);
        buildMarkovModel(backgroundModel, bgSequences);
        countKmers(kmerCounts1, sequence1, score.kmerSize);
        countKmers(kmerCounts2, sequence2, score.kmerSize);

        unsigned nvals = length(kmerCounts1);  // Number of kmers
        int len1 = 0;
        int len2 = 0;

        for (unsigned l = 0; l < nvals; l++)
        {
            len1 += kmerCounts1[l];
            len2 += kmerCounts2[l];
        }
        String<TValue> probabilities;
        resize(probabilities, nvals, missing);
        for (unsigned i = 0; i < nvals; ++i)
        {
            TValue p_w = 1.0;  // Probability of kmer
            TValue variance = 0.0;
            String<TUnmaskedAlphabet> w;
            unhash(w, i, score.kmerSize);
            p_w = emittedProbability(backgroundModel, w);
            variance = ((TValue) pow(((TValue) len1 * len2), 0.5)) * p_w;
            TValue variance1 = 0.0;
            TValue variance2 = 0.0;

            variance1 = pow(len1 * p_w, 0.5);
            variance2 = pow(len2 * p_w, 0.5);

            // Calculate standardised kmer Count
            if ((variance > pow(10.0, -10)) && (variance < pow(10.0, 10)))
            {
                if (p_w > 0)
                {
                    TValue stCount1 = (kmerCounts1[i] - p_w * len1) / variance1;
                    TValue stCount2 = (kmerCounts2[i] - p_w * len2) / variance2;
                    result += stCount1 * stCount2;
                 }
            }
        }
    }
}

}  // namespace seqan

#endif  // SEQAN_INCLUDE_SEQAN_ALIGNMENT_FREE_AF_D2STAR_ORIGINAL_H_
