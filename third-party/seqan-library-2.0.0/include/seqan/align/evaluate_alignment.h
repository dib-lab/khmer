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
// Compute alignment score given a pairwise alignment.
// ==========================================================================

#ifndef INCLUDE_SEQAN_ALIGN_EVALUATE_ALIGNMENT_H_
#define INCLUDE_SEQAN_ALIGN_EVALUATE_ALIGNMENT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class AlignmentStats
// ----------------------------------------------------------------------------

/*!
 * @class AlignmentStats
 * @headerfile <seqan/align.h>
 * @brief Statistics about a tabular alignment.
 *
 * @signature struct AlignmentStats;
 *
 * @see computeAlignmentStats
 *
 * @fn AlignmentStats::AlignmentStats
 * @brief Constructor
 *
 * @signature AlignmentStats::AlignmentStats();
 *
 * All members are initialized to <tt>0</tt>.
 *
 * @var unsigned AlignmentStats::numGapOpens;
 * @brief Number of gap open events.
 *
 * @var unsigned AlignmentStats::numGapExtensions;
 * @brief Number of gap extension events.
 *
 * @var unsigned AlignmentStats::numInsertions;
 * @brief Number of gaps in reference relative to query.
 *
 * @var unsigned AlignmentStats::numDeletions;
 * @brief Number of gaps in query relative to reference.
 *
 * @var unsigned AlignmentStats::numMatches;
 * @brief Number of match (identity) events.
 *
 * @var unsigned AlignmentStats::numMismatches;
 * @brief Number of mismatch (not identity) events.
 *
 * @var unsigned AlignmentStats::numPositiveScores;
 * @brief Number of residues aligned with positive score (0 is counted as positive).
 *
 * @var unsigned AlignmentStats::numNegativeScores;
 * @brief Number of residues aligned with negative score.
 *
 * @var float AlignmentStats::alignmentSimilarity;
 * @brief The resulting alignment percent similarity (positive).
 *
 * @var float AlignmentStats::alignmentIdentity;
 * @brief The resulting alignment percent identity (match).
 *
 * @var int AlignmentStats::alignmentScore;
 * @brief The resulting alignment score.
 */

struct AlignmentStats
{
    // Number of gap opens/gap extensions.
    unsigned numGapOpens;
    unsigned numGapExtensions;
    // Number of insertions and deletions.
    unsigned numInsertions;
    unsigned numDeletions;
    // Number of matches, mismatches.
    unsigned numMatches;
    unsigned numMismatches;
    // Number of aligned residues with positive/negative scores.
    unsigned numPositiveScores;
    unsigned numNegativeScores;

    // the alignment identity and similarity scores
    float alignmentSimilarity;
    float alignmentIdentity;

    // The alignment score.
    int alignmentScore;

    AlignmentStats() : numGapOpens(0), numGapExtensions(0), numInsertions(0), numDeletions(0),
                       numMatches(0), numMismatches(0), numPositiveScores(0), numNegativeScores(0),
                       alignmentSimilarity(0.0), alignmentIdentity(0.0), alignmentScore(0)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

/*!
 * @fn AlignmentStats#clear
 * @brief Resets all members to <tt>0</tt>.
 *
 * @signature void clear(stats);
 *
 * @param[in,out] stats AlignmentStats object to clear.
 */

inline
void clear(AlignmentStats & stats)
{
    stats.numGapOpens = 0;
    stats.numGapExtensions = 0;
    stats.numInsertions = 0;
    stats.numDeletions = 0;
    stats.numMatches = 0;
    stats.numMismatches = 0;
    stats.numPositiveScores = 0;
    stats.numNegativeScores = 0;
    stats.alignmentSimilarity = 0.0;
    stats.alignmentIdentity = 0.0;
    stats.alignmentScore = 0;
}

// ----------------------------------------------------------------------------
// Function computeAlignmentStats()
// ----------------------------------------------------------------------------

/*!
 * @fn computeAlignmentStats
 * @headerfile <seqan/align.h>
 * @brief Compute alignment statistics.
 *
 * @signature TScoreVal computeAlignmentStats([stats, ]align, scoringScheme);
 *
 * @param[out] stats The @link AlignmentStats @endlink object to store alignment statistics in.
 * @param[in]  align The @link Align @endlink object to score.
 * @param[in]  score The @link Score @endlink object to use for the scoring scheme.
 *
 * @return TScoreVal The score value of the alignment, of the same type as the value type of <tt>scoringScheme</tt>
 *
 * @see AlignmentStats
 *
 * @section Examples
 *
 * @include demos/align/compute_alignment_stats.cpp
 *
 * The output is as follows:
 *
 * @include demos/align/compute_alignment_stats.cpp.stdout
 */

template <typename TSource, typename TAlignSpec, typename TScoreVal, typename TScoreSpec>
TScoreVal computeAlignmentStats(AlignmentStats & stats,
                                Align<TSource, TAlignSpec> const & align,
                                Score<TScoreVal, TScoreSpec> const & scoringScheme)
{
    SEQAN_ASSERT_EQ_MSG(length(rows(align)), 2u, "Only works with pairwise alignments.");
    SEQAN_ASSERT_EQ_MSG(length(row(align, 0)), length(row(align, 1)), "Invalid alignment!");
    clear(stats);

    typedef Align<TSource, TAlignSpec> const TAlign;
    typedef typename Row<TAlign>::Type TGaps;
    typedef typename Iterator<TGaps, Standard>::Type TGapsIter;
    typedef typename Value<typename Source<TGaps>::Type>::Type TAlphabet;

    // Get iterators.
    TGapsIter it0 = begin(row(align, 0));
    TGapsIter itEnd0 = end(row(align, 0));
    TGapsIter it1 = begin(row(align, 1));
    TGapsIter itEnd1 = end(row(align, 1));

    // State whether we have already opened a gap.
    bool isGapOpen0 = false, isGapOpen1 = false;

    for (; it0 != itEnd0 && it1 != itEnd1; ++it0, ++it1)
    {
        if (isGap(it0))
        {
            if (isGapOpen0)
            {
                stats.numGapOpens += 1;
                stats.alignmentScore += scoreGapOpen(scoringScheme);
            }
            else
            {
                stats.numGapExtensions += 1;
                stats.alignmentScore += scoreGapExtend(scoringScheme);
            }
            stats.numDeletions += 1;
            isGapOpen0 = true;
        }
        else
        {
            isGapOpen0 = false;
        }

        if (isGap(it1))
        {
            if (!isGapOpen1)
            {
                stats.numGapOpens += 1;
                stats.alignmentScore += scoreGapOpen(scoringScheme);
            }
            else
            {
                stats.numGapExtensions += 1;
                stats.alignmentScore += scoreGapExtend(scoringScheme);
            }
            stats.numInsertions += 1;
            isGapOpen1 = true;
        }
        else
        {
            isGapOpen1 = false;
        }

        if (!isGap(it0) && !isGap(it1))
        {
            // Compute the alignment score and register in stats.
            TAlphabet c0 = *it0, c1 = *it1;
            TScoreVal scoreVal = score(scoringScheme, c0, c1);
            stats.alignmentScore += scoreVal;
            // Register other statistics.
            bool isMatch = (c0 == c1);
            bool isPositive = (scoreVal >= 0);
            stats.numMatches += isMatch;
            stats.numMismatches += !isMatch;
            stats.numPositiveScores += isPositive;
            stats.numNegativeScores += !isPositive;
        }
    }
    SEQAN_ASSERT(it0 == itEnd0);
    SEQAN_ASSERT(it1 == itEnd1);

    // Finally, compute the alignment similarity from the various counts
    float alignmentLength = static_cast<float>(length(row(align, 0)));
    stats.alignmentSimilarity = 100.0 * static_cast<float>(stats.numPositiveScores) / alignmentLength;
    stats.alignmentIdentity = 100.0 * static_cast<float>(stats.numMatches) / alignmentLength;

    return stats.alignmentScore;
}

template <typename TGaps, typename TAlignSpec, typename TScoreVal, typename TScoreSpec>
TScoreVal computeAlignmentStats(Align<TGaps, TAlignSpec> const & align,
                                Score<TScoreVal, TScoreSpec> const & scoringScheme)
{
    AlignmentStats stats;
    (void)stats;
    return computeAlignmentStats(stats, align, scoringScheme);
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_EVALUATE_ALIGNMENT_H_
