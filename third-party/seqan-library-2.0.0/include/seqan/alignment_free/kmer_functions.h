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
// This file contains helper functions to count words in sequences and to
// calculate probabilities and variances of word occurrences.
// ==========================================================================

// TODO (goeke) const could be added below for the input variables but the function value() in matrix_base (align) is not defined for const.  Similarly, the function emittedProbabilty is not defined for const in statistics_markov_model.h

#ifndef SEQAN_INCLUDE_SEQAN_ALIGNMENT_FREE_KMER_FUNCTIONS_H_
#define SEQAN_INCLUDE_SEQAN_ALIGNMENT_FREE_KMER_FUNCTIONS_H_

namespace seqan {

template <typename TAlphabet>
struct UnmaskedAlphabet_
{
    typedef TAlphabet Type;
};

template <>
struct UnmaskedAlphabet_<Dna5>
{
    typedef Dna Type;
};

template <typename TAlphabet>
struct UnmaskedAlphabet_<const TAlphabet>
{
    typedef const typename UnmaskedAlphabet_<TAlphabet>::Type Type;
};

/*!
 * @fn countKmers
 * @headerfile <seqan/alignment_free.h>
 * @brief Counts kmers in a sequence. Optionally, a background model is returned.
 *
 * @signature void countKmers(kmerCounts, sequence, k);
 * @signature void countKmers(kmerCounts, bgFrequencies, sequence, k);
 * @signature void countKmers(kmerCounts, bgModel, sequence, k);
 *
 * @param[out] kmerCounts    @link String @endlink of <tt>unsigned</tt> with kmer counts for every k-mer.
 * @param[out] bgFrequencies @link String @endlink of background frequencies (<tt>double</tt>) representing the model.
 * @param[out] bgModel       @link MarkovModel @endlink to use.
 * @param[in]  sequence      @link String @endlink (sequence) where k-mers are counted.
 * @param[in]  k             k-mer length (<tt>unsigned</tt>).
 *
 * k-mers overlapping masked (aka 'N') letters are not counted in case of Dna5Strings.  A Bernoulli or Markov Model
 * can be choosen as a background model.
 *
 * @section Examples
 *
 * Calculate the alignment free sequence similarity o two masked DNA sequences.
 *
 * @code{.cpp}
 * using namespace seqan;
 * // Masked sequence, we do not want to count words overlapping 'N'
 * Dna5String sequenceDna5 =
 *     "TAGGTTTTCCGAAAAGGTAGCAACTTTACGTGATCAAACCTCTGACGGGGTTTTCCCCGTCGAAATTGGGTG"
 *     "TTTCTTGTCTTGTTCTCACTTGGGGCATCTCCGTCAAGCCAAGAAAGTGCTCCCTGGATTCTGTTGCTAACG"
 *     "AGTCTCCTCTGCATTCCTGCTTGACTGATTGGGCGGACGGGGTGTCCACCTGACGCTGAGTATCGCCGTCAC"
 *     "GGTGCCACATGTCTTATCTATTCAGGGATCAGAATTTATTCAGGAAATCAGGAGATGCTACACTTGGGTTAT"
 *     "CGAAGCTCCTTCCAAGGCGTAGCAAGGGCGACTGAGCGCGTAAGCTCTAGATCTCCTCGTGTTGCAACTACA"
 *     "CGCGCGGGTCACTCGAAACACATAGTATGAACTTAACGACTGCTCGTACTGAACAATGCTGAGGCAGAAGAT"
 *     "CGCAGACCAGGCATCCCACTGCTTGAAAAAACTATNNNNCTACCCGCCTTTTTATTATCTCATCAGATCAAG";
 *
 * String<unsigned> kmerCounts;
 * unsigned k = 2;  // Count all 2-mers
 * countKmers(kmerCounts, sequenceDna5, k);
 *
 * for(unsigned i = 0; i<16; ++i)       // Print the 2-mer counts
 *     std::cout<<kmerCounts[i]<<"\n";  // 34 times AA; 30 times AC; 28 times AG; ...
 *
 *
 * String<double> nucleotideFrequencies;  // Defines a Bernoulli model for DNA sequences.
 * // Count all 2-mers and save the nucleotide frequencies
 * countKmers(kmerCounts, nucleotideFrequencies, sequenceDna5, k);
 *
 * for(unsigned i = 0; i<4; ++i)          // Print the nucleotide frequencies
 *     std::cout << nucleotideFrequencies[i] << "\n";
 * // => p(A) = 0.238; p(C) = 0.254; p(G) = 0.238; p(T) = 0.27;
 *
 * MarkovModel<Dna, double>  backgroundModel(1);  // Markov model of order 1
 * // Count all 2-mers and return a Markov model
 * countKmers(kmerCounts, backgroundModel, sequenceDna5, k);
 * std::cout<<backgroundModel.transition;  // Print the transition matrix of the markov model
 * @endcode
 *
 * @see alignmentFreeComparison
 * @see calculateProbability
 * @see calculateVariance
 * @see calculateCovariance
 * @see stringToStringSet
 * @see MarkovModel
 */

/*
 * Function to count kmers, Ns are not considered
 */
template <typename TString>
void countKmers(String<unsigned> & kmerCounts, TString const & sequence, unsigned const k)
{
    typedef typename Value<TString>::Type               TAlphabet;
    typedef typename UnmaskedAlphabet_<TAlphabet>::Type TUnmaskedAlphabet;
    typedef typename Iterator<TString const>::Type      TIterator;
    typedef typename Position<TIterator>::Type          TPosition;
    typedef Shape<TUnmaskedAlphabet, SimpleShape>       TShape;
    // Declare variables
    TShape myShape;  // Shape, length can be changed (kmer_length)
    resize(myShape, k);
    // Calculate the number of kmers, length of count vector
    int kmerNumber = _intPow((unsigned)ValueSize<TUnmaskedAlphabet>::VALUE, weight(myShape));
    clear(kmerCounts);
    resize(kmerCounts, kmerNumber, 0);
    TIterator itSequence = begin(sequence);
    int counterN = 0;

    // Check for any N that destroys the first kmers
    unsigned j = k - 1;
    for (TPosition i = position(itSequence); i <= j; ++i)
    {
        if (_repeatMaskValue(sequence[i]))
        {
            counterN = i + 1;
        }
    }
    for (; itSequence <= (end(sequence) - k); ++itSequence)
    {
        // Check if there is a "N" at the end of the new kmer
        if (_repeatMaskValue(value(itSequence + (k - 1))))
            counterN = k;  // Do not consider any kmer covering this "N"

        // If there is no "N" overlapping with the current kmer, count it
        if (counterN <= 0)
        {
            unsigned hashValue = hash(myShape, itSequence);
            ++kmerCounts[hashValue];
        }
        counterN--;
    }
}

/*
 * Function to count kmers and background nucleotide frequencies, Ns are not considered
 * (for zero order background model)
 */
template <typename TValueBG, typename TString>
void countKmers(String<unsigned> & kmerCounts, String<TValueBG> & backgroundFrequencies, TString const & sequence, unsigned const k)
{
    typedef typename Value<TString>::Type                    TAlphabet;
    typedef typename UnmaskedAlphabet_<TAlphabet>::Type      TUnmaskedAlphabet;
    typedef typename Iterator<TString const>::Type           TIterator;
    typedef typename Iterator<String<TValueBG> >::Type       TIteratorTStringBG;
    typedef typename Position<TIterator>::Type               TPosition;
    typedef Shape<TUnmaskedAlphabet, SimpleShape>            TShape;
    unsigned alphabetSize = ValueSize<TUnmaskedAlphabet>::VALUE;

    // Declare variables
    TShape myShape;         // Shape, length can be changed (kmer_length)
    TShape myShapeBG;       // Shape for background, set to markovlen+1, here zero order only
    resize(myShape, k);
    resize(myShapeBG, 1);   // Markov model of zero order (count background frequencies)

    // Calculate number of kmers/ length of count vector, respectively background vector
    unsigned kmerNumber = _intPow(alphabetSize, k);
    unsigned kmerNumberBG = alphabetSize;  // Zero order model for DNA sequences (Bernoulli model)
    clear(kmerCounts);
    resize(kmerCounts, kmerNumber, 0);
    resize(backgroundFrequencies, kmerNumberBG, (TValueBG) 0);
    TIterator itSequence = begin(sequence);

    int counterN = 0;  // Counter that counts how many kmers are effected by a N
    int counterNbg = 0;   // Counter for background model (different shape size)

    // Check for any N that destroys the first kmers
    unsigned j = k - 1;
    for (TPosition i = position(itSequence); i <= j; ++i)
    {
        if (_repeatMaskValue(sequence[i]))
            counterN = i + 1;
    }

    int sumBG = 0;  // Count the number of nucleotides for the nucleotide frequency calculation (Ns are not considered anymore).
    for (; itSequence <= (end(sequence) - k); ++itSequence)
    {
        // Check if there is a "N" at the end of the new kmer
        if (_repeatMaskValue(value(itSequence + (k - 1))))
        {
            counterN = k;  // Do not consider any kmer covering this "N"
        }
        // If there is no "N" overlapping with the current kmer, count it.
        if (counterN <= 0)
        {
            unsigned hashValue = hash(myShape, itSequence);
            ++kmerCounts[hashValue];
        }
        // Check if there is a "N" at the end of the new background word, here single letters only.
        if (_repeatMaskValue(value(itSequence)))
        {
            counterNbg = 1;
        }
        if (counterNbg <= 0)
        {
            unsigned hashValueBG = hash(myShapeBG, itSequence);
            backgroundFrequencies[hashValueBG] += 1.0;
            ++sumBG;
        }
        counterN--;
        counterNbg--;
    }
    // The background counts are updated until the last base is covered.
    for (; itSequence < end(sequence); ++itSequence)
    {
        if (_repeatMaskValue(value(itSequence)))
        {
            counterNbg = 1;
        }
        if (counterNbg <= 0)
        {
            unsigned hashValueBG = hash(myShapeBG, itSequence);
            ++backgroundFrequencies[hashValueBG];
            ++sumBG;
        }
        counterNbg--;
    }
    // Normalise the background counts to obtain the nucleotide frequencies (Bernoulli model of DNA sequences).
    TIteratorTStringBG itBackground = begin(backgroundFrequencies);
    for (; itBackground < end(backgroundFrequencies); ++itBackground)
        if (sumBG != 0)
            value(itBackground) /= ((TValueBG) sumBG);
}

/*
 * Function to count kmers and build a background markov model with masked sequences
 */
template <typename TString, typename TAlphabetBG, typename TValue>
void countKmers(String<unsigned> & kmerCounts, MarkovModel<TAlphabetBG, TValue> & backgroundModel, TString const & sequence, unsigned k)
{
    //typedef typename Value<TString>::Type                   TAlphabet;
    //typedef typename UnmaskedAlphabet_<TAlphabet>::Type     TUnmaskedAlphabet;
    typedef typename Iterator<TString const, Rooted>::Type  TIterator;
    //typedef typename Iterator<String<int>, Rooted>::Type    TIteratorInt;
    typedef typename Position<TIterator>::Type              TPosition;
    typedef Shape<TAlphabetBG, SimpleShape>                 TShape;

    // Declare variables
    TShape myShape;  // Shape, length can be changed (kmer_length)
    resize(myShape, k);
    // Only consider kmers without N
    int kmerNumber = _intPow((unsigned)ValueSize<TAlphabetBG>::VALUE, weight(myShape));
    clear(kmerCounts);
    resize(kmerCounts, kmerNumber, 0);

    // Create sequence set for the markov model, if Ns occur, the sequence is split and Ns are removed
    StringSet<String<TAlphabetBG> > seqSetMM;

    TIterator itSeq = begin(sequence);

    // Check for any N that destroys the first kmers
    unsigned j = (k - 1);
    for (TPosition i = position(itSeq); i <= j; ++i)
    {
        if (_repeatMaskValue(sequence[i]))
        {
            if ((i - position(itSeq)) > 0)
            {
                appendValue(seqSetMM, infix(sequence, position(itSeq), i));
            }
            goFurther(itSeq, i + 1 - position(itSeq));
            j = i + k - 1;
        }
    }

    int counterN = 0;
    TPosition startSplitSequence = position(itSeq);  // The position of possible start of a sequence after NNs is stored to split sequences.

    for (; itSeq <= (end(sequence) - k); ++itSeq)
    {
        if (_repeatMaskValue(value(itSeq + (k - 1))))
        {
            counterN = k;

            if (((position(itSeq) + k - 1) > startSplitSequence))
                appendValue(seqSetMM, infix(sequence, startSplitSequence, (position(itSeq) + k - 1)));

            startSplitSequence = (position(itSeq) + k);  // Position after N, possible start
        }
        if (counterN <= 0)
        {
            unsigned hashValue = hash(myShape, itSeq);
            ++kmerCounts[hashValue];
        }

        counterN--;
    }
    // Create a stringSet, needed to create the Markov model
    if ((position(itSeq) + k - 1) > startSplitSequence)
    {
        appendValue(seqSetMM, infix(sequence, startSplitSequence, (position(itSeq) + k - 1)));
    }
    // Build background Markov model
    buildMarkovModel(backgroundModel, seqSetMM);
}

/*!
 * @fn calculateProbability
 * @headerfile <seqan/alignment_free.h>
 * @brief Calculates the probability of a sequence given a Bernoulli model.
 *
 * @signature void calculateProbability(probability, sequence, bgFrequencies);
 *
 * @param[out] probability  Probability (<tt>double</tt>) of the sequence given the model.
 * @param[in] sequence      @link String @endlink, usually of Dna characters.
 * @param[in] bgFrequencies @link String @endlink of background frequencies (<tt>double</tt>) representing the model.
 *
 * @section Examples
 *
 * Calculate the probability for the word CCCAAGTTT with <i>p(A) = p(T) = 0.3</i> and <i>p(C) = p(G) = 0.2</i>.
 *
 * @code{.cpp}
 * using namespace seqan;
 * double p = 0.0;
 * DnaString word = "CCCAAGTTT";
 * String<double> model;
 * resize(model, 4);
 * model[0] = 0.3;  // p(A)
 * model[1] = 0.2;  // p(C)
 * model[2] = 0.2;  // p(G)
 * model[3] = 0.3;  // p(T)
 * calculateProbability(p, word, model);  // p = 3.888e-06
 * @endcode
 *
 * @see calculateVariance
 * @see alignmentFreeComparison
 * @see calculateCovariance
 * @see countKmers
 */

template <typename TValue, typename TString, typename TStringBG>
void calculateProbability(TValue & probability, TString const & sequence, TStringBG const & backgroundFrequencies)
{
    typedef typename Iterator<TString const, Rooted>::Type  TIteratorTString;

    TIteratorTString itSequence = begin(sequence);
    probability = (TValue) 1;
    for (; itSequence < end(sequence); ++itSequence)
        probability *= backgroundFrequencies[ordValue(*itSequence)];
}

/*!
 * @fn calculateVariance
 * @headerfile <seqan/alignment_free.h>
 * @brief Calculates the variance for the number of word occurrences of a word in a sequence of length n given a
 *        background model.
 *
 * @signature void calculateVariance(variance, word, bgFrequencies, n);
 * @signature void calculateVariance(variance, word, bgModel, n);
 *
 * @param[out] variance     Variance of the number of occurrences of the word in a sequence of length n given the
 *                          model; <tt>double</tt>.
 * @param[in] word          @link String @endlink, usually of Dna to compute variance for.
 * @param[in] bgFrequencies @link String @endlink of bg frequencies representing the model.
 * @param[in] bgModel       @link MarkovModel @endlink to use.
 * @param[in] n             Length of the sequence where the occurrences of word are counted, <tt>int</tt>.
 *
 * Calculates the variance for the number of word occurrences of a word in a sequence of length n given a background
 * model (Markov model or Bernoulli model). The formula is obtained from (Robin et al., 2005).
 *
 * @section References
 *
 * Robin, S., Rodolphe, F., and Schbath, S.  (2005). DNA, Words and Models. Cambridge University Press. See Jonathan
 * Goeke et al (to appear) for details on the implementation.
 *
 * @section Examples
 *
 * Calculate the variance for the number of occurrences of CAAGTC in a sequence of length 10000bp with
 * <i>p(A) = p(T) = 0.3</i> and <i>p(C) = p(G) = 0.2</i>.
 *
 * @code{.cpp}
 * using namespace seqan;
 * double var = 0.0;
 * int n = 10000;
 * DnaString word = "CAAGTC";
 * String<double> model;
 * resize(model, 4);
 * model[0] = 0.3;  // p(A)
 * model[1] = 0.2;  // p(C)
 * model[2] = 0.2;  // p(G)
 * model[3] = 0.3;  // p(T)
 * calculateVariance(var, word, model, n);  // var = 2.16
 * @endcode
 *
 * Estimate a Markov model on a set of sequences and calculate the variance for the number of occurrences of the word
 * CAAGTC in a sequence of length 10000bp.
 *
 * @code{.cpp}
 * using namespace seqan;
 * double var = 0.0;
 * int n = 10000;
 * DnaString word = "CAAGTC";
 * StringSet<DnaString> sequences;
 * appendValue(sequences, "CAGAAAAAAACACTGATTAACAGGAATAAGCAGTTTACTTATTTTGGGCCTGGGACCCGTGTCTCTAATTTAATTAGGTGATCCCTGCGAAGTTTCTCCA");
 * MarkovModel<Dna, double> model(0);  // Bernoulli model
 * model.build(sequences);
 * calculateVariance(var, word, model, n);  // var = 2.16
 * MarkovModel<Dna, double> model1(1);  // First order Markov model
 * model1.build(sequences);
 * calculateVariance(var, word, model1, n);  // var = 1.69716
 * @endcode
 *
 * @see calculateProbability
 * @see calculateCovariance
 * @see MarkovModel
 * @see alignmentFreeComparison
 * @see calculatePeriodicity
 * @see countKmers
 * @see calculateOverlapIndicator
 */

template <typename TValue, typename TString, typename TStringBG>
void calculateVariance(TValue & variance, TString const & word, TStringBG const & backgroundFrequencies, int const n)
{
    typedef typename Value<TString>::Type                       TAlphabet;
    typedef typename Value<TStringBG>::Type                     TValueBG;
    typedef typename Iterator<String<int>, Rooted>::Type        TIteratorInt;

    int l = length(word);
    TValueBG p_w;
    calculateProbability(p_w, word, backgroundFrequencies);

    String<int> periodicity;
    calculatePeriodicity(periodicity, word, word);
    variance = (TValue) (n - l + 1) * p_w;
    for (TIteratorInt i = begin(periodicity); i < end(periodicity); ++i)
    {
        TValueBG p_clump;
        TValueBG p_tmp;
        calculateProbability(p_tmp, word, backgroundFrequencies);
        String<TAlphabet> wordPrefix = prefix(word, value(i));
        calculateProbability(p_clump, wordPrefix, backgroundFrequencies);
        p_clump *= p_tmp;
        variance += (TValue) 2 * (n - l + 1 - value(i)) * p_clump;
    }
    variance += (TValue) p_w * p_w * (n - 2 * n * l + 3 * l * l - 4 * l + 1);
}


template <typename TValue, typename TSpec, typename TAlphabet>
void calculateVariance(TValue & variance, String<TAlphabet, TSpec> const & word, MarkovModel<TAlphabet, TValue> /*const*/ & bgModel, int const n)
{
    typedef typename Iterator<String<int>, Rooted>::Type        TIteratorInt;

    int l = length(word);
    TValue p_w;
    p_w = emittedProbability(bgModel, word);
    String<int> periodicity;
    calculatePeriodicity(periodicity, word, word);
    variance = (TValue) (n - l + 1) * p_w;
    for (TIteratorInt i = begin(periodicity); i < end(periodicity); ++i)
    {
        TValue p_clump;
        String<TAlphabet> clump = prefix(word, value(i));
        append(clump, word);
        p_clump = emittedProbability(bgModel, clump);
        variance += (TValue) 2 * (n - l + 1 - value(i)) * p_clump;
    }
    variance += (TValue) p_w * p_w * (n - 2 * n * l + 3 * l * l - 4 * l + 1);
}

/*!
 * @fn calculateCovariance
 * @headerfile <seqan/alignment_free.h>
 * @brief Calculates the covariance for the number of word occurrences for two words in a sequence of length n, given a
 *        background model.
 *
 * @signature void calculateCovariance(covariance, word1, word2, bgFrequencies, n);
 * @signature void calculateCovariance(covariance, word1, word2, bgModel, n);
 *
 * @param[out] covariance    Variance of the number of occurrences of the word in a sequence of length n given the
 *                           model, <tt>double</tt>.
 * @param[in] word1          @link String @endlink, usually of Dna.
 * @param[in] word2          @link String @endlink, usually of Dna.
 * @param[in] bgFrequencies  @link String @endlink of <tt>double</tt> with the background frequencies representing
 * @param[in] bgModel        @link MarkovModel @endlink to use.
 * @param[in] n              Length of the sequence where the occurrences of word are counted, <tt>int</tt>.
 *
 * Calculates the covariance for the number of word occurrences for two words in a sequence of length n given a
 * background model (Markov model or Bernoulli model). The covariance is influenced by the property of words to overlap,
 * for example, the words ATAT and TATA have a high covariance since they are likely to overlap. The formula is based on
 * (Robin et al., 2005).
 *
 * @section References
 *
 * Robin, S., Rodolphe, F., and Schbath, S.  (2005). DNA, Words and Models. Cambridge University Press. See Jonathan
 * Goeke et al (to appear) for details on the implementation.
 *
 * @section Examples
 *
 * Calculate the covariance for the number of occurrences of ATATAT and TATATA in a sequence of length 10000bp with
 * <i>p(A) = p(T) = 0.3</i> and <i>p(C) = p(G) = 0.2</i>.
 *
 * @code{.cpp}
 * using namespace seqan;
 * double covar = 0.0;
 * int n = 10000;
 * DnaString word1 = "ATATAT";
 * DnaString word2 = "TATATA";
 * String<double> model;
 * resize(model, 4);
 * model[0] = 0.3;  // p(A)
 * model[1] = 0.2;  // p(C)
 * model[2] = 0.2;  // p(G)
 * model[3] = 0.3;  // p(T)
 * calculateCovariance(covar, word1, word2, model, n);  // covar = 4.74
 * @endcode
 *
 * Estimate a Markov model on a set of sequences and calculate the covariance for the number of occurrences of ATATAT
 * and TATATA in a sequence of length 10000bp.
 *
 * @code{.cpp}
 * using namespace seqan;
 * double covar = 0.0;
 * int n = 10000;
 * DnaString word1 = "ATATAT";
 * DnaString word2 = "TATATA";
 * StringSet<DnaString> sequences;
 * appendValue(sequences, "CAGCACTGATTAACAGGAATAAGCAGTTTACTTCTGTCAGAATATTGGGCATATATA"
 *                        "CTGGGACCCGTGTAATACTCTAATTTAATTAGGTGATCCCTGCGAAGTCTCCA");
 * MarkovModel<Dna, double> modelMM0(0);  // Bernoulli model
 * modelMM0.build(sequences);
 * calculateCovariance(covar, word1, word2, modelMM0, n);  // covar = 4.74
 * MarkovModel<Dna, double> modelMM1(1);  // First order Markov model
 * modelMM1.build(sequences);
 * calculateCovariance(covar, word1, word2, modelMM1, n);  // covar = 13.1541
 * @endcode
 *
 * @see calculateProbability
 * @see calculateVariance
 * @see MarkovModel
 * @see alignmentFreeComparison
 * @see calculatePeriodicity
 * @see countKmers
 * @see calculateOverlapIndicator
 */

template <typename TValue, typename TString, typename TStringBG>
void calculateCovariance(TValue & covariance, TString const & word1, TString const & word2, TStringBG const & backgroundFrequencies, int const n)
{
    if (word1 == word2)
    {
        calculateVariance(covariance, word1, backgroundFrequencies, n);
        return;
    }
    typedef typename Value<TString>::Type                   TAlphabet;
    typedef typename Value<TStringBG>::Type                 TValueBG;
    typedef typename Iterator<String<int>, Rooted>::Type    TIteratorInt;

    covariance = 0;
    int l1 = length(word1);
    TValueBG p_w1;
    calculateProbability(p_w1, word1, backgroundFrequencies);
    String<int> periodicity1;
    calculatePeriodicity(periodicity1, word1, word2);
    for (TIteratorInt i = begin(periodicity1); i < end(periodicity1); ++i)
    {
        TValueBG p_clump;
        TValueBG p_tmp;
        calculateProbability(p_tmp, word2, backgroundFrequencies);
        String<TAlphabet> wordPrefix = prefix(word1, value(i));
        calculateProbability(p_clump, wordPrefix, backgroundFrequencies);
        p_clump *= p_tmp;
        covariance += (TValue) (n - l1 + 1 - value(i)) * p_clump;
    }

    int l2 = length(word2);
    TValueBG p_w2;
    calculateProbability(p_w2, word2, backgroundFrequencies);
    String<int> periodicity2;
    calculatePeriodicity(periodicity2, word2, word1);
    for (TIteratorInt i = begin(periodicity2); i < end(periodicity2); ++i)
    {
        TValueBG p_clump;
        TValueBG p_tmp;
        calculateProbability(p_tmp, word1, backgroundFrequencies);
        String<TAlphabet> wordPrefix = prefix(word2, value(i));
        calculateProbability(p_clump, wordPrefix, backgroundFrequencies);
        p_clump *= p_tmp;
        covariance += (TValue) (n - l2 + 1 - value(i)) * p_clump;
    }
    covariance += (TValue) p_w1 * p_w2 * (n - 2 * n * l1 + 3 * l1 * l1 - 4 * l1 + 1);
}

template <typename TValue, typename TSpec, typename TAlphabet>
void calculateCovariance(TValue & covariance, String<TAlphabet, TSpec> const & word1, String<TAlphabet, TSpec> const & word2, MarkovModel<TAlphabet, TValue> /*const*/ & bgModel, int const n)
{
    if (word1 == word2)
    {
        calculateVariance(covariance, word1, bgModel, n);
        return;
    }
    typedef typename Iterator<String<int>, Rooted>::Type TIteratorInt;

    covariance = 0;
    int l1 = length(word1);
    TValue p_w1;
    p_w1 = emittedProbability(bgModel, word1);
    String<int> periodicity1;
    calculatePeriodicity(periodicity1, word1, word2);  // word2 is right
    for (TIteratorInt i = begin(periodicity1); i < end(periodicity1); ++i)
    {
        TValue p_clump;
        String<TAlphabet> clump = prefix(word1, value(i));
        append(clump, word2);

        p_clump = emittedProbability(bgModel, clump);

        covariance += (TValue) (n - l1 + 1 - value(i)) * p_clump;
    }
    TValue p_w2;
    p_w2 = emittedProbability(bgModel, word2);
    String<int> periodicity2;
    calculatePeriodicity(periodicity2, word2, word1);
    for (TIteratorInt i = begin(periodicity2); i < end(periodicity2); ++i)
    {
        TValue p_clump;
        String<TAlphabet> clump = prefix(word2, value(i));
        append(clump, word1);

        p_clump = emittedProbability(bgModel, clump);

        covariance += (TValue) (n - l1 + 1 - value(i)) * p_clump;
    }
    covariance += (TValue) p_w1 * p_w2 * (n - 2 * n * l1 + 3 * l1 * l1 - 4 * l1 + 1);
}

/*!
 * @fn calculatePeriodicity
 * @headerfile <seqan/alignment_free.h>
 * @brief Calculate word periodicity (indicator for overlaps)
 *
 * @signature void calculatePeriodicity(periodicity, word1, word2);
 *
 * @param[out] periodicity String of <tt>int</tt> values giving the periodicity (overlap indicator) of
 *                         <tt>word1</tt> and <tt>word2</tt>.
 * @param[int] word1       String, usually of Dna characters.
 * @param[int] word2       String, usually of Dna characters.
 *
 * Calculate word periodicity (indicator for overlaps) for two words.
 *
 * @section Examples
 *
 * Calculate the periodicity of two words (At which positions can they overlap?)
 *
 * @code{.cpp}
 * using namespace seqan;
 * DnaString word1 = "ATATA";
 * DnaString word2 = "TATAT";
 * String<int> periodicity;
 * calculatePeriodicity(periodicity, word1, word2);
 * for(unsigned i = 0; i < length(periodicity); ++i)  // Print the periodicity
 *     std::cout << periodicity[i] << "\t";
 *
 * // periodocity[0] = 1:
 * // i =     01234
 * // word1 = ATATA
 * // word2 = -TATAT
 *
 * // periodocity[1] = 3:
 * // i =     01234
 * // word1 = ATATA
 * // word2 = ---TATAT
 * @endcode
 *
 * @see calculateVariance
 * @see calculateCovariance
 * @see calculateOverlapIndicator
 * @see alignmentFreeComparison
 */

template <typename TString>
void calculatePeriodicity(String<int> & periodicity, TString const & word1, TString const & word2)
{
    typedef typename Value<TString>::Type                   TAlphabet;
    //typedef typename Iterator<TString const, Rooted>::Type  TIterator;
    typedef typename Size<TString>::Type                    TSize;

    TSize length1 = length(word1);
    TSize length2 = length(word2);
    for (TSize i = 1; i < length1; ++i)
    {
        String<TAlphabet> my_suffix = suffix(word1, i);  // Overlap of suffix of word1 with prefix of word2
        TSize my_min = std::min(length2, (length1 - i));
        String<TAlphabet> my_prefix = prefix(word2, my_min);
        if (my_suffix == my_prefix)
        {
            appendValue(periodicity, i);
        }
    }
}

/*!
 * @fn calculateOverlapIndicator
 * @headerfile <seqan/alignment_free.h>
 * @brief Calculate word overlaps: <tt>epsilon(word1, word2) = 1</tt> where <tt>word2[j] = word1[j+p] for
 *        all j = 1..(k-p)</tt>.
 *
 * @signature void calculateOverlapIndicator(epsilon, word1, word2);
 *
 * @param[out] epsilon String of int giving the periodicity (overlap indicator) of word1 and word2.
 * @param[in]  word1   String (for example a DNA sequence).
 * @param[in]  word2   String (for example a DNA sequence).
 *
 * Calculate the indicator for overlaps of two words. The formula is based on (Robin et al., 2005)
 *
 * @section References
 *
 * Robin, S., Rodolphe, F., and Schbath, S. (2005). DNA, Words and Models.  Cambridge University Press. See Jonathan
 * Goeke et al (to appear) for details on the implementation.
 *
 * @section Examples
 *
 * Calculate the overlap indicator (epsilon) for two words
 *
 * @code{.cpp}
 * using namespace seqan;
 * DnaString word1 = "ATATA";
 * DnaString word2 = "TATAT";
 * String<int> epsilon;
 * calculateOverlapIndicator(epsilon, word1, word2);
 * for(unsigned i = 0; i < length(epsilon); ++i)
 *     std::cout << epsilon[i] << "\t";
 * // epsilon =         01010:
 * // word1             ATATA
 * // word2 overlap 1:  -TATAT
 * // word2 overlap 2:  ---TATAT
 * @endcode
 *
 * @see calculateVariance
 * @see calculateCovariance
 * @see calculatePeriodicity
 * @see alignmentFreeComparison
 */

template <typename TString>
void calculateOverlapIndicator(String<int> & epsilon, TString const & word1, TString const & word2)
{
    typedef typename Value<TString>::Type                   TAlphabet;
    //typedef typename Iterator<TString const, Rooted>::Type  TIterator;
    typedef typename Size<TString>::Type                    TSize;

    TSize length1 = length(word1);
    TSize length2 = length(word2);
    clear(epsilon);
    resize(epsilon, length1, 0);
    for (TSize i = 0; i < length1; ++i)
    {
        String<TAlphabet> my_suffix = suffix(word1, length1 - i - 1);  // Overlap of suffix of word1 with prefix of word2
        TSize my_min = std::min(length2, i + 1);
        String<TAlphabet> my_prefix = prefix(word2, my_min);
        if (my_suffix == my_prefix)
            epsilon[i] = 1;
    }
}

/*!
 * @fn stringToStringSet
 * @headerfile <seqan/alignment_free.h>
 * @brief Transform a String into a StringSet containing this String.
 *
 * @signature void stringToStringSet(stringSet, string);
 * @signature void stringToStringSet(dnaStringSet, dna5String);
 *
 * @param[out] stringSet    @link StringSet @endlink to create with one sequence.
 * @param[in]  string       @link String @endlink to create the string set of.
 * @param[out] dnaStringSet @link StringSet @endlink of @link String Strings @endlink over the alphabet @link Dna @endlink.
 * @param[in]  dna5String   @link String @endlink over the alphabet @link Dna5 @endlink to convert.
 *
 * @note The second variant removes all N characters from the @link Dna5String @endlink.
 *
 * @section Examples
 *
 * Transform a masked DNA sequence into a set of sequences with all masked parts removed.
 *
 * @code{.cpp}
 * using namespace seqan;
 * Dna5String sequenceDna5 =
 *     "NNNNNNTTTCCGAAAAGGTANNNNNGCAACTTTANNNCGTGATCAAAGTTTTCCCCGTCGAAATTGGGNNTG";
 * StringSet<DnaString> sequencesDna;
 * stringToStringSet(sequencesDna, sequenceDna5);
 * // Print the masked sequence
 * std::cout<<sequenceDna5<<"\n";
 * // Print the sequence with the masked parts removed
 * for(unsigned i = 0; i < length(sequencesDna); ++i)
 *     std::cout<<sequencesDna[i]<<"\n";
 * // sequencesDna[0] = "TTTCCGAAAAGTA"
 * // sequencesDna[1] = "GCAACTTTA"
 * // sequencesDna[2] = "CGTGATCAAAGTTTTCCCCGTCGAAATTGGG"
 * // sequencesDna[3] = "TG"
 * @endcode
 *
 * @see MarkovModel
 * @see alignmentFreeComparison
 * @see cutNs
 * @see countKmers
 */

template <typename TString>
void
stringToStringSet(StringSet<TString> & stringSet, TString const & sequence)
{
    resize(stringSet, 1);
    stringSet[0] = sequence;
}

inline void
stringToStringSet(StringSet<String<Dna> > & dnaStringSet, String<Dna5> const & sequence)
{
    typedef Iterator<String<Dna5> const, Rooted>::Type  TIterator;
    typedef Position<TIterator>::Type                   TPosition;

    TIterator itSeq = begin(sequence);
    // Check for any N that destroys the first kmers
    unsigned j = 0;
    for (TPosition i = position(itSeq); i <= j; ++i)
    {
        if (sequence[i] == 'N')
        {
            if ((i - position(itSeq)) > 0)
                appendValue(dnaStringSet, infix(sequence, position(itSeq), i));
            goFurther(itSeq, i + 1 - position(itSeq));
            j = i;
        }
    }
    int counterN = 0;
    TPosition startSplitSequence = position(itSeq);  // The position of possible starts of a sequence after Ns is stored to split the sequence.
    for (; itSeq <= (end(sequence) - 1); ++itSeq)
    {
        if (value(itSeq) == 'N')
        {
            counterN = 1;
            if (((position(itSeq)) > startSplitSequence))
            {
                appendValue(dnaStringSet, infix(sequence, startSplitSequence, position(itSeq)));
            }
            startSplitSequence = (position(itSeq) + 1);  // Position after N, possible start
        }
        counterN--;
    }
    // Create the stringSet, the stringSet can be used to create a Markov model
    if (position(itSeq) > startSplitSequence)
        appendValue(dnaStringSet, infix(sequence, startSplitSequence, position(itSeq)));
}

/*!
 * @fn cutNs
 * @headerfile <seqan/alignment_free.h>
 * @brief Cut out all masked sequences from a Dna5String.
 *
 * @signature void cutNs(sequenceCut, sequence);
 *
 * @param[out] sequenceCut Dna5String similar to sequence with all Ns cut out.
 * @param[in]  sequence    Masked DNA sequence.
 *
 * This function concatenates the nonmasked parts of the sequence, thereby changing the word content. If you want to
 * remove the masked parts of a sequence without concatenation, use stringToStringSet.
 *
 * @section Examples
 *
 * Transform a masked DNA sequence into an unmasked sequences with all masked parts cut out
 *
 * @code{.cpp}
 * using namespace seqan;
 * Dna5String sequenceMasked =
 *     "NNNNNNTTTCCGAAAAGGTANNNNNGCAACTTTANNNCGTGATCAAAGTTTTCCCCGTCGAAATTGGGNNTG";
 * Dna5String sequenceMaskedPartsRemoved;
 * cutNs(sequenceMaskedPartsRemoved, sequenceMasked);
 * // Print the masked sequence
 * std::cout<<sequenceMasked<<"\n";
 * // Print the sequence with the masked parts removed
 * std::cout<<sequenceMaskedPartsRemoved<<"\n";
 * // sequenceMasked =
 * // "NNNNNNTTTCCGAAAAGGTANNNNNGCAACTTTANNNCGTGATCAAAGTTTTCCCCGTCGAAATTGGGNNTG"
 * // sequenceMaskedPartsRemoved =
 * // "TTTCCGAAAAGGTAGCAACTTTACGTGATCAAAGTTTTCCCCGTCGAAATTGGGTG"
 * @endcode
 *
 * @see MarkovModel
 * @see stringToStringSet
 * @see alignmentFreeComparison
 */

inline void
cutNs(String<Dna5> & sequenceCut, String<Dna5> const & sequence)
{
    typedef Iterator<String<Dna5> const, Rooted>::Type  TIterator;
    typedef Position<TIterator>::Type                   TPosition;

    sequenceCut = "";
    TIterator itSeq = begin(sequence);

    // Check for any N that destroys the first kmers
    unsigned j = 0;
    for (TPosition i = position(itSeq); i <= j; ++i)
    {
        if (sequence[i] == 'N')
        {
            if ((i - position(itSeq)) > 0)
                sequenceCut += infix(sequence, position(itSeq), i);
            goFurther(itSeq, i + 1 - position(itSeq));
            j = i;
        }
    }
    int counterN = 0;
    TPosition startSplitSequence = position(itSeq);  // The position of possible starts of a sequence after Ns is stored to split the sequence.
    for (; itSeq <= (end(sequence) - 1); ++itSeq)
    {
        if (value(itSeq) == 'N')
        {
            counterN = 1;
            if (((position(itSeq)) > startSplitSequence))
                sequenceCut += infix(sequence, startSplitSequence, position(itSeq));
            startSplitSequence = (position(itSeq) + 1);  // Position after N, possible start
        }
        counterN--;
    }
    // Create the sequence with any N cut out.
    if (position(itSeq) > startSplitSequence)
    {
        sequenceCut += infix(sequence, startSplitSequence, position(itSeq));
    }
}

}  // namespace seqan

#endif  // SEQAN_INCLUDE_SEQAN_ALIGNMENT_FREE_KMER_FUNCTIONS_H_
