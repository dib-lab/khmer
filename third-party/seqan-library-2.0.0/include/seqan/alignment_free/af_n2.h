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
// This header contains the implementation of the N2 score for alignment free
// sequence comparison with word neighbourhood counts
//
// See: Goeke et al, to appear.
//
// These functions can be called with alignmentFreeComparison().
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGNMENT_FREE_AF_N2_H_
#define SEQAN_INCLUDE_SEQAN_ALIGNMENT_FREE_AF_N2_H_

namespace seqan {

/*
 * This function returns a string with indices indicating which k-mer is the
 * reverse complement k-mer: i=revComIndex[revComIndex[i]]
 */
inline void _initialiseRevComIndex(String<unsigned> & revComIndex, unsigned const k)
{
    unsigned myLength = (unsigned)pow(4.0, (int)k);
    resize(revComIndex, myLength, 0);
    Shape<Dna, SimpleShape> myShape;
    resize(myShape, k);
    for (unsigned i = 0; i < myLength; ++i)
    {
        String<Dna> w;
        unhash(w, i, k);
        DnaStringReverseComplement wRC(w);
        unsigned hashValue = hash(myShape, begin(wRC));
        revComIndex[i] = hashValue;
    }

}

/*
 * This function returns a stringSet with strings of indices indicating which
 * k-mers belong to the word neighbourhood for every k-mer (all k-mers with
 * one mismatch)
 */
inline void _initialiseKmerNeighbourhood(StringSet<String<unsigned> > & kmerNeighbourhood,
                                  unsigned const k, bool const revCom,
                                  String<unsigned> const & revComIndex)
{
    unsigned myLength = (unsigned)pow(4.0, (int)k);
    Shape<Dna, SimpleShape> myShape;
    resize(myShape, k);
    resize(kmerNeighbourhood, myLength);
    for (unsigned i = 0; i < myLength; ++i)
    {
        resize(kmerNeighbourhood[i], 1, i);

        String<Dna> w;
        unhash(w, i, k);
        if ((revComIndex[i] != i) && (revCom == true))
        {
            appendValue(kmerNeighbourhood[i], revComIndex[i]);
        }
        for (unsigned j = 0; j < k; ++j)
        {
            for (unsigned l = 0; l < 4; ++l)
            {
                String<Dna> wTMP;
                wTMP = w;
                if (wTMP[j] != l)
                {
                    wTMP[j] = l;
                    unsigned hashValue = hash(myShape, begin(wTMP));
                    // Check for double word occurrences
                    bool duplicate = false;
                    if (revCom == true)
                    {

                        for (unsigned n = 0; n < length(kmerNeighbourhood[i]); ++n)
                        {
                            if ((hashValue) == kmerNeighbourhood[i][n])
                            {
                                duplicate = true;
                                break;
                            }
                        }

                    }

                    if (duplicate == false)
                    {
                        appendValue(kmerNeighbourhood[i], hashValue);
                        if (revCom == true)
                        {
                            if (revComIndex[hashValue] != hashValue)
                            {
                                appendValue(kmerNeighbourhood[i], revComIndex[hashValue]);
                            }
                        }
                    }
                }
            }
        }
    }
}

/*
 * _alignmentFreeComparison is called by alignmentFreeComparison() (see alignment_free_comparison.h)
 */
template <typename TValue, typename TStringSet>
void _alignmentFreeComparison(Matrix<TValue, 2> & scoreMatrix,
                              TStringSet const & sequenceSet,
                              AFScore<N2> const & score)
{



    typedef typename Value<TStringSet>::Type                            TString;
    typedef typename Value<TString>::Type                               TAlphabet;
    typedef typename UnmaskedAlphabet_<TAlphabet>::Type                 TUnmaskedAlphabet;
    typedef typename Iterator<TStringSet const>::Type                   TIteratorSet;
    typedef typename Iterator<StringSet<String<double> > >::Type        TIteratorSetDouble;


    // Initialise the reverse complement hash table
    String<unsigned> revComIndex;
    StringSet<String<unsigned> > kmerNeighbourhood;
    _initialiseRevComIndex(revComIndex, score.kmerSize);
    if (score.revCom == "both_strands")
    {
        _initialiseKmerNeighbourhood(kmerNeighbourhood, score.kmerSize, true, revComIndex);
    }
    else
    {
        _initialiseKmerNeighbourhood(kmerNeighbourhood, score.kmerSize, false, revComIndex);
    }

    unsigned seqNumber = length(sequenceSet);

    setLength(scoreMatrix, 0, seqNumber);
    setLength(scoreMatrix, 1, seqNumber);
    resize(scoreMatrix, (TValue) 0);

    StringSet<String<double> > standardisedKmerCounts;
    resize(standardisedKmerCounts, seqNumber);
    // Count all kmers and all background nucleotide frequencies and store them in StringSets
    TIteratorSetDouble itStandardisedKmerCounts = begin(standardisedKmerCounts);
    TIteratorSet itSeqSet = begin(sequenceSet);
    for (; itSeqSet < end(sequenceSet); ++itSeqSet)
    {

        _standardiseCounts(*itStandardisedKmerCounts, revComIndex, kmerNeighbourhood, *itSeqSet, score);
        if(score.verbose)
        {
            std::cout << "\n" << position(itSeqSet);
        }
        ++itStandardisedKmerCounts;
    }

    if (score.norm == true) // Normalise the score so that sequence-self-comparisons are always 1
    {
        itStandardisedKmerCounts = begin(standardisedKmerCounts);
        for (; itStandardisedKmerCounts < end(standardisedKmerCounts); ++itStandardisedKmerCounts)
        {
            TValue normValue = 0.0;
            for (unsigned i = 0; i < length(value(itStandardisedKmerCounts)); ++i)
            {
                normValue += value(itStandardisedKmerCounts)[i] * value(itStandardisedKmerCounts)[i];
            }
            for (unsigned i = 0; i < length(value(itStandardisedKmerCounts)); ++i)
            {
                value(itStandardisedKmerCounts)[i] /= sqrt(normValue);
            }

        }
    }

    if (!(score.outputFile == ""))
    {
        std::ofstream myfile;
        myfile.open(toCString(score.outputFile));
        for (unsigned i = 0; i < length(standardisedKmerCounts[0]); ++i)
        {
            String<TUnmaskedAlphabet> w;
            unhash(w, i, score.kmerSize);
            myfile << "\t" << w;
        }
        myfile << "\n";
        for (unsigned seqIndex = 0; seqIndex < seqNumber; ++seqIndex)
        {
            myfile << "Seq" << seqIndex;
            for (unsigned i = 0; i < length(standardisedKmerCounts[seqIndex]); ++i)
            {
                myfile << "\t" << standardisedKmerCounts[seqIndex][i];
            }
            myfile << "\n";
        }
        myfile.close();
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
            _alignmentFreeCompareCounts(value(scoreMatrix, rowIndex, colIndex), revComIndex, standardisedKmerCounts[rowIndex], standardisedKmerCounts[colIndex], score);
            value(scoreMatrix, colIndex, rowIndex) = value(scoreMatrix, rowIndex, colIndex);  // Copy symmetric entries
        }
    }
}

/*
 * Calculate pairwise score given the counts of all kmers
 */
template <typename TValue, typename TString>
void
_alignmentFreeCompareCounts(TValue & result,
                            String<unsigned> const revComIndex,
                            TString const & kmerCounts1,
                            TString const & kmerCounts2,
                            AFScore<N2> const & score)
{
    typedef typename Iterator<TString const, Rooted>::Type    TIteratorTString;

    TIteratorTString it1 = begin(kmerCounts1);
    TIteratorTString it2 = begin(kmerCounts2);
    result = 0.0;
    TValue resultRC = 0.0;
    for (; it1 < end(kmerCounts1); ++it1)
    {
        result += (TValue)(value(it1) * value(it2));
        // Computation of the reverse complement strand score
        if ((score.revCom != "") && (score.revCom != "both_strands"))
        {
            unsigned hashValue = revComIndex[position(it1)];
            resultRC += (TValue)(value(it1) * kmerCounts2[hashValue]);
        }
        ++it2;
    }

    if (score.revCom == "mean")
    {
        result = (TValue) (resultRC + result) / 2;
    }
    else if (score.revCom == "max")
    {
        result = std::max(resultRC, result);
    }
    else if (score.revCom == "min")
    {
        result = std::min(resultRC, result);
    }
}

/*
 * count kmers and standardise count vectors for Dna5 and markov model background
 */
template <typename TString, typename TSequence>
void _standardiseCounts(TString & standardisedCounts,
                        String<unsigned> const & revComIndex,
                        StringSet<String<unsigned> > const & kmerNeighbourhood,
                        TSequence const & sequence,
                        AFScore<N2> const & score)
{
    typedef typename Value<TSequence>::Type                     TAlphabet;
    typedef typename UnmaskedAlphabet_<TAlphabet>::Type         TUnmaskedAlphabet;
    typedef typename Value<TString>::Type                       TValue;
    typedef typename Iterator<String<unsigned>, Rooted>::Type   TIteratorUnsigned;
    typedef typename Iterator<TString, Rooted>::Type            TIteratorTString;

    unsigned alphabetSize = ValueSize<TUnmaskedAlphabet>::VALUE;

    // Save all word covariances which are computed in covariance Matrix to avoid double computations
    Matrix<TValue, 2> covarianceMatrix;
    TValue missing = -pow(10.0, 10);
    if (score.mismatches > 0)
    {
        setLength(covarianceMatrix, 0, pow((double)alphabetSize, (int)score.kmerSize));
        setLength(covarianceMatrix, 1, pow((double)alphabetSize, (int)score.kmerSize));
        resize(covarianceMatrix, missing);
    }

    // Note that there is some code below that looks like copy-and-paste.  However, pulling this out into another
    // function is the only way to get rid of the duplicate lines since we use different types.  After some discussion,
    // weese, goeke and holtgrew agreed that it is probably easier to read and maintain this way than to spread the code
    // over to one more function.
    if (score.bgModelOrder == 0)
    {
        // ----------------------------------------------------------------------
        // Order 0 Background Model
        // ----------------------------------------------------------------------

        String<unsigned> kmerCounts;
        String<double> backgroundFrequencies;
        countKmers(kmerCounts, backgroundFrequencies, sequence, score.kmerSize);
        int nvals = length(kmerCounts);  // Number of kmers
        int len1 = 0;
        for (int l = 0; l < nvals; l++)
        {
            len1 += kmerCounts[l];
        }
        resize(standardisedCounts, nvals, (TValue) 0.0);

        // String of TValue to store the word probabilites p_w
        String<TValue> probabilities;
        resize(probabilities, nvals, missing);

        TIteratorUnsigned itCounts;
        TIteratorTString itStandardisedCounts;

        itCounts = begin(kmerCounts);
        itStandardisedCounts = begin(standardisedCounts);

        for (; itCounts < end(kmerCounts); ++itCounts)
        {
            // Temporary counter for mismatch kmer counting
            TValue counterTMP = 0;
            TValue p_w = 1;  // Probability of kmer

            String<TUnmaskedAlphabet> w;
            unhash(w, (unsigned)position(itCounts), score.kmerSize);
            calculateProbability(p_w, w, backgroundFrequencies);
            TValue variance = 0;
            if ((score.mismatches == 1))  // Mismatch  score calculation
            {
                p_w = 0;
                // The first word in the kmerNeighbourhood is the kmer itself, it is weighted normally
                // Sum of all entries in the covariance matrix. only once computed
                unsigned wordHash = position(itCounts);
                unsigned wordRCHash = revComIndex[wordHash];

                for (unsigned row = 0; row < length(kmerNeighbourhood[wordHash]); ++row)
                {
                    unsigned wordRowHash = kmerNeighbourhood[wordHash][row];
                    // The kmer itself is weighted normally
                    if (wordRowHash == wordHash)  // The first word in the kmerNeighbourhood is the kmer itself, it is weighted normally
                    {
                        counterTMP += (TValue) kmerCounts[wordRowHash];
                    }
                    else if ((score.revCom == "both_strands") && (wordRowHash == wordRCHash))
                    {
                        counterTMP += ((TValue) kmerCounts[wordRowHash]);
                    }
                    else
                    {
                        counterTMP += ((TValue) kmerCounts[wordRowHash]) * score.mismatchWeight;
                    }
                    String<Dna> wMM1;
                    unhash(wMM1, wordRowHash, score.kmerSize);

                    for (unsigned col = row; col < length(kmerNeighbourhood[wordHash]); ++col)
                    {
                        unsigned wordColHash = kmerNeighbourhood[wordHash][col];
                        if (value(covarianceMatrix, wordColHash, wordRowHash) == missing)
                        {
                            String<Dna> wMM2;

                            unhash(wMM2, wordColHash, score.kmerSize);
                            calculateCovariance(value(covarianceMatrix, wordColHash, wordRowHash), wMM1, wMM2, backgroundFrequencies, (len1 + score.kmerSize - 1));
                            value(covarianceMatrix, wordRowHash, wordColHash) = value(covarianceMatrix, wordColHash, wordRowHash);
                        }
                        if (row == col)  // Variance of weighted variables
                        {
                            if ((wordRowHash == wordHash) || (score.revCom == "both_strands" && (wordRowHash == wordRCHash)))  // The variance of the kmer is counted full
                            {
                                variance += value(covarianceMatrix, wordRowHash, wordColHash);
                            }
                            else
                            {
                                variance += pow(score.mismatchWeight, 2) * value(covarianceMatrix, wordRowHash, wordColHash);  // Calculate weighted variances
                            }
                        }
                        // The covariance of the kmer and the reverse complement is weighted full
                        else if ((score.revCom == "both_strands") && (((wordRowHash == wordHash) && (wordColHash == wordRCHash)) || ((wordRowHash == wordRCHash) && (wordColHash == wordHash))))
                        {
                            variance += (2.0) * value(covarianceMatrix, wordRowHash, wordColHash);
                        }
                        else if ((wordRowHash == wordHash || wordColHash == wordHash) || (score.revCom == "both_strands" && (wordRowHash == wordRCHash || wordColHash == wordRCHash)))  // The covariance is weighted half
                        {
                            variance += (2.0) * score.mismatchWeight * value(covarianceMatrix, wordRowHash, wordColHash);
                        }
                        else  // The covariance is weighted^2
                        {
                            variance += (2.0) * pow(score.mismatchWeight, 2) * value(covarianceMatrix, wordRowHash, wordColHash);
                        }
                    }
                    if (probabilities[wordRowHash] == missing)
                    {
                        calculateProbability(probabilities[wordRowHash], wMM1, backgroundFrequencies);
                    }
                    if (wordRowHash == wordHash)  //Weight the probabilities and expected values, normal weight for the kmer itself
                    {
                        p_w += probabilities[wordRowHash];
                    }
                    else if ((score.revCom == "both_strands") && (wordRowHash == wordRCHash))  // Weight the probabiliets and expected values, normal weight for the reverse complement kmer itself
                    {
                        p_w += probabilities[wordRowHash];
                    }
                    else
                    {
                        p_w += score.mismatchWeight * probabilities[wordRowHash];
                    }
                }
                variance = pow(variance, 0.5);
            }  // End of mismatch calculation
            else if (score.revCom == "both_strands")
            {
                TValue variance1;
                TValue variance2;
                TValue covariance;
                String<Dna> wRC;
                unhash(wRC, (unsigned)  revComIndex[(unsigned)position(itCounts)], score.kmerSize);
                calculateVariance(variance1, w, backgroundFrequencies, (len1 + score.kmerSize - 1));
                calculateVariance(variance2, wRC, backgroundFrequencies, (len1 + score.kmerSize - 1));
                calculateCovariance(covariance, w, wRC, backgroundFrequencies, (len1 + score.kmerSize - 1));
                variance = pow((variance1 + variance2 + (2.0) * covariance), 0.5);
                TValue p_wRC = 1;  // Probability of the reverse complement kmer
                calculateProbability(p_wRC, wRC, backgroundFrequencies);
                p_w += p_wRC;
            }
            else
            {
                calculateVariance(variance, w, backgroundFrequencies, (len1 + score.kmerSize - 1));
                variance = pow(variance, 0.5);
            }
            if ((variance > pow(10.0, -10)) && (variance < pow(10.0, 10)))
            {
                if (p_w > 0)
                {
                    if (score.mismatches > 0)
                    {
                        value(itStandardisedCounts) = ((TValue) ((TValue) counterTMP) - p_w * ((TValue)len1)) / variance;
                    }
                    else if (score.revCom == "both_strands")
                    {
                         value(itStandardisedCounts) = ((TValue) ((TValue) value(itCounts) + kmerCounts[revComIndex[(unsigned)position(itCounts)]]) - p_w * ((TValue)len1)) / variance;
                    }
                    else
                    {
                        value(itStandardisedCounts) = ((TValue) ((TValue) value(itCounts)) - p_w * ((TValue)len1)) / variance;
                    }
                }
            }
            ++itStandardisedCounts;
        }
    }
    else
    {
        // ----------------------------------------------------------------------
        // Higher Order Background Model
        // ----------------------------------------------------------------------

        String<unsigned> kmerCounts;
        MarkovModel<TUnmaskedAlphabet, TValue> backgroundModel(score.bgModelOrder);
        countKmers(kmerCounts, backgroundModel, sequence, score.kmerSize);

        int nvals = length(kmerCounts);  // Number of kmers
        int len1 = 0;
        for (int l = 0; l < nvals; l++)
        {
            len1 += kmerCounts[l];
        }
        resize(standardisedCounts, nvals, (TValue) 0.0);
        String<TValue> probabilities;
        resize(probabilities, nvals, missing);
        TIteratorUnsigned itCounts;
        TIteratorTString itStandardisedCounts;
        itCounts = begin(kmerCounts);
        itStandardisedCounts = begin(standardisedCounts);

        for (; itCounts < end(kmerCounts); ++itCounts)
        {
            TValue p_w = 1;  // Probability of kmer
            TValue variance = 0;
            String<TUnmaskedAlphabet> w;
            unhash(w, (unsigned)position(itCounts), score.kmerSize);
            p_w = emittedProbability(backgroundModel, w);

            TValue counterTMP = 0.0;
            if ((score.mismatches == 1))  // Start of mismatch calculations
            {
                p_w = 0;
                // The first word in the kmerNeighbourhood is the kmer itself, it is weighted normally
                // Sum of all entries in the covariance matrix, computed and stored dynamically
                unsigned wordHash = position(itCounts);
                unsigned wordRCHash = revComIndex[wordHash];

                for (unsigned row = 0; row < length(kmerNeighbourhood[wordHash]); ++row)
                {
                    unsigned wordRowHash = kmerNeighbourhood[wordHash][row];
                    // The kmer itself is weighted normally
                    if (wordRowHash == wordHash)  // The first word in the kmerNeighbourhood is the kmer itself, it is weighted normally
                    {
                        counterTMP += (TValue) kmerCounts[wordRowHash];
                    }
                    else if ((score.revCom == "both_strands") && (wordRowHash == wordRCHash))
                    {
                        counterTMP += ((TValue) kmerCounts[wordRowHash]);
                    }
                    else
                    {
                        counterTMP += ((TValue) kmerCounts[wordRowHash]) * score.mismatchWeight;
                    }
                    String<Dna> wMM1;
                    unhash(wMM1, wordRowHash, score.kmerSize);
                    for (unsigned col = row; col < length(kmerNeighbourhood[wordHash]); ++col)
                    {
                        unsigned wordColHash = kmerNeighbourhood[wordHash][col];
                        if (value(covarianceMatrix, wordColHash, wordRowHash) == missing)
                        {
                            String<Dna> wMM2;
                            unhash(wMM2, wordColHash, score.kmerSize);
                            calculateCovariance(value(covarianceMatrix, wordColHash, wordRowHash), wMM1, wMM2, backgroundModel, (len1 + score.kmerSize - 1));
                            value(covarianceMatrix, wordRowHash, wordColHash) = value(covarianceMatrix, wordColHash, wordRowHash);
                        }
                        if (row == col)  // Variance of weighted variables
                        {
                            if ((wordRowHash == wordHash) || (score.revCom == "both_strands" && (wordRowHash == wordRCHash)))  // The variance of the kmer is counted full
                            {
                                variance += value(covarianceMatrix, wordRowHash, wordColHash);
                            }
                            else
                            {
                                variance += pow(score.mismatchWeight, 2) * value(covarianceMatrix, wordRowHash, wordColHash);
                            }
                        }
                        // The covariance of the kmer and the reverse complement is weighted full
                        else if ((score.revCom == "both_strands") && (((wordRowHash == wordHash) && (wordColHash == wordRCHash)) || ((wordRowHash == wordRCHash) && (wordColHash == wordHash))))
                        {
                            variance += (2.0) * value(covarianceMatrix, wordRowHash, wordColHash);
                        }
                        else if ((wordRowHash == wordHash || wordColHash == wordHash) || (score.revCom == "both_strands" && (wordRowHash == wordRCHash || wordColHash == wordRCHash)))  // The covariance is weighted half
                        {
                            variance += (2.0) * score.mismatchWeight * value(covarianceMatrix, wordRowHash, wordColHash);
                        }
                        else  // The covariance is weighted^2
                        {
                            variance += (2.0) * pow(score.mismatchWeight, 2) * value(covarianceMatrix, wordRowHash, wordColHash);
                        }
                    }
                    if (probabilities[wordRowHash] == missing)
                    {
                        probabilities[wordRowHash] = emittedProbability(backgroundModel, wMM1);
                    }
                    if (wordRowHash == wordHash) //Weight the probabiliets and expected values, normal weight for the kmer itself
                    {
                        p_w += probabilities[wordRowHash];
                    }
                    else if ((score.revCom == "both_strands") && (wordRowHash == wordRCHash)) // Weight the probabiliets and expected values, normal weight for the reverse complement kmer itself
                    {
                        p_w += probabilities[wordRowHash];
                    }
                    else
                    {
                        p_w += score.mismatchWeight * probabilities[wordRowHash];
                    }
                }
                variance = pow(variance, 0.5);  // Calculate the standard deviation
            }  // End of mismatch calculations
            else if (score.revCom == "both_strands")
            {
                TValue variance1;
                TValue variance2;
                TValue covariance;
                String<Dna> wRC;
                unhash(wRC, (unsigned)revComIndex[(unsigned)position(itCounts)], score.kmerSize);
                calculateVariance(variance1, w, backgroundModel, (len1 + score.kmerSize - 1));
                calculateVariance(variance2, wRC, backgroundModel, (len1 + score.kmerSize - 1));
                calculateCovariance(covariance, w, wRC, backgroundModel, (len1 + score.kmerSize - 1));
                variance = pow((variance1 + variance2 + (2.0) * covariance), 0.5);
                TValue p_wRC = 1;   // Probability of the reverse complement kmer
                p_wRC = emittedProbability(backgroundModel, wRC);
                p_w += p_wRC;
            }
            else
            {
                calculateVariance(variance, w, backgroundModel, (len1 + score.kmerSize - 1));
                variance = pow(variance, 0.5);
            }
            if ((variance > pow(10.0, -10)) && (variance < pow(10.0, 10)))
            {
                if (p_w > 0)
                {
                    if (score.mismatches > 0)
                    {
                        value(itStandardisedCounts) = ((TValue) ((TValue) counterTMP) - p_w * ((TValue)len1)) / variance;
                    }
                    else if (score.revCom == "both_strands")
                    {
                        value(itStandardisedCounts) = ((TValue) ((TValue) value(itCounts) + kmerCounts[revComIndex[(unsigned)position(itCounts)]]) - p_w * ((TValue)len1)) / variance;
                    }
                    else
                    {
                        value(itStandardisedCounts) = ((TValue) ((TValue) value(itCounts)) - p_w * ((TValue)len1)) / ((TValue) variance);
                    }
                }
                ++itStandardisedCounts;
            }
        }
    }
}

}  // namespace seqan

#endif  // SEQAN_INCLUDE_SEQAN_ALIGNMENT_FREE_AF_N2_H_
