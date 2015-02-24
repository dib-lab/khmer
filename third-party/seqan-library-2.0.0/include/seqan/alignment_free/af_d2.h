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
// This header contains the implementation of the D2 score for alignment free
// sequence comparison (inner product of kmer counts).
//
// These functions can be called with alignmentFreeComparison().
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGNMENT_FREE_AF_D2_H_
#define SEQAN_INCLUDE_SEQAN_ALIGNMENT_FREE_AF_D2_H_

namespace seqan {

/*
 * _alignmentFreeComparison is called by alignmentFreeComparison() (see alignment_free_comparison.h)
 */
template <typename TStringSet, typename TValue>
void _alignmentFreeComparison(Matrix<TValue, 2> & scoreMatrix,
                              TStringSet const & sequenceSet,
                              AFScore<D2> const & score)
{
    //typedef typename Value<TStringSet>::Type                        TString;
    //typedef typename Value<TString>::Type                           TAlphabet;
    typedef typename Iterator<TStringSet const>::Type               TIteratorSet;
    typedef typename Iterator<StringSet<String<unsigned> > >::Type  TIteratorSetInt;
    //typedef Matrix<TValue, 2>                                       TMatrix;

    unsigned seqNumber = length(sequenceSet);

    // Resize the scoreMatrix
    setLength(scoreMatrix, 0, seqNumber);
    setLength(scoreMatrix, 1, seqNumber);
    resize(scoreMatrix, (TValue) 0);

    StringSet<String<unsigned> > kmerCounts;
    resize(kmerCounts, seqNumber);

    // Count all kmers
    TIteratorSetInt itKmerCounts = begin(kmerCounts);
    TIteratorSet itSeqSet = begin(sequenceSet);

    for (; itSeqSet < end(sequenceSet); ++itSeqSet)
    {
        countKmers(value(itKmerCounts), value(itSeqSet), score.kmerSize);
        ++itKmerCounts;
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
            _alignmentFreeCompareCounts(value(scoreMatrix, rowIndex, colIndex), kmerCounts[rowIndex], kmerCounts[colIndex], score);
            value(scoreMatrix, colIndex, rowIndex) = value(scoreMatrix, rowIndex, colIndex);  // Copy symmetric entries
        }
    }

}

/*
 * Calculate pairwise score given the counts of all kmers
 */
template <typename TValue>
void
_alignmentFreeCompareCounts(TValue & result,
                            String<unsigned> const & kmerCounts1,
                            String<unsigned> const & kmerCounts2,
                            AFScore<D2> const & /*score*/)
{
    typedef typename Iterator<String<unsigned> const>::Type TIteratorInt;

    TIteratorInt it1 = begin(kmerCounts1);
    TIteratorInt it2 = begin(kmerCounts2);

    result = 0;
    for (; it1 < end(kmerCounts1); ++it1)
    {
        result += value(it1) * value(it2);
        ++it2;
    }
}

}  // namespace seqan

#endif  // SEQAN_INCLUDE_SEQAN_ALIGNMENT_FREE_AF_D2_H_
