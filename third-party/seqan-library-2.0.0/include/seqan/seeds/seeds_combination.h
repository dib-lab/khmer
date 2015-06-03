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
// Algorithms for combining (i.e. merging and chaining) seeds.
// ==========================================================================

// TODO(holtgrew): All the Nothing()'s should not be part of the public interface.

#ifndef SEQAN_SEEDS_SEEDS_COMBINATION_H_
#define SEQAN_SEEDS_SEEDS_COMBINATION_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

/*!
 * @defgroup LocalChainingTags
 * @brief Tags for selecting local chaining algorithms.
 *
 * @tag LocalChainingTags#Merge
 * @headerfile <seqan/seeds.h>
 * @brief Try to merge with existing seed.
 *
 * @signature typedef Tag<Merge_> Merge;
 *
 * @tag LocalChainingTags#Chaos
 * @headerfile <seqan/seeds.h>
 * @brief Chain to other seed using CHAOS chaining condition.
 *
 * @signature typedef Tag<Chaos_> Chaos;
 *
 * @tag LocalChainingTags#SimpleChain
 * @headerfile <seqan/seeds.h>
 * @brief Chain to other seed by simple chaining.
 *
 * @signature typedef Tag<SimpleChain_> SimpleChain;
 *
 * @tag LocalChainingTags#Single
 * @headerfile <seqan/seeds.h>
 * @brief Force adding without chaining.
 *
 * @signature typedef Tag<Single_> Single;
 */

// TODO(holtgrew): Stream-line tags to Merge, ChaosChain, SimpleChain?
struct Merge_;
typedef Tag<Merge_> Merge;

struct Chaos_;
typedef Tag<Chaos_> Chaos;

struct SimpleChain_;
typedef Tag<SimpleChain_> SimpleChain;

struct Single_;
typedef Tag<Single_> Single;

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

// Returns true iff b can be merged into a where a is the one to the
// upper left, b the one to the lower right.
template <typename TSeedSpec, typename TSeedConfig, typename TThreshold, typename TBandwidth>
inline bool
_seedsCombineable(Seed<TSeedSpec, TSeedConfig> const & seedA,
                  Seed<TSeedSpec, TSeedConfig> const & seedB,
                  TThreshold const & maxDiagonalDistance,
                  TBandwidth const & /*maxBandwidth*/,
                  Merge const &)
{
    // TODO(holtgrew): TThreshold could be Position<TSeed>::Type.
    // b has to be right of a for the two seeds to be mergeable.
    if (beginPositionH(seedB) < beginPositionH(seedA) || beginPositionV(seedB) < beginPositionV(seedA))
        return false;
    // If the two seeds do not overlap, they cannot be merged.
    if (beginPositionH(seedB) > endPositionH(seedA) || beginPositionV(seedB) > endPositionV(seedA))
        return false;
    // If the distance between the diagonals exceeds the threshold
    // then the seeds cannot be merged.
    typedef typename MakeUnsigned_<TThreshold>::Type TUnsignedThreshold;
    if (static_cast<TUnsignedThreshold>(_abs(endDiagonal(seedA) - beginDiagonal(seedB))) > static_cast<TUnsignedThreshold>(maxDiagonalDistance))
        return false;
    // Otherwise, the seeds can be merged.
    return true;
}


// Returns true iff b can be simple-chained to a where a is the one to
// the upper left, b the one to the lower right.
template <typename TSeedSpec, typename TSeedConfig, typename TThreshold, typename TBandwidth>
inline bool
_seedsCombineable(Seed<TSeedSpec, TSeedConfig> const & seedA,
                  Seed<TSeedSpec, TSeedConfig> const & seedB,
                  TThreshold const & maxGapSize,
                  TBandwidth const & /*maxBandwidth*/,
                  SimpleChain const &)
{
    // TODO(holtgrew): We should be able to configure whether we want to have Manhattan, euclidean, minimal edit distance, for seeds.
    // TODO(holtgrew): TThreshold could be Position<TSeed>::Type.

    // b has to be right of a for the two seeds to be chainable.
    if (beginPositionH(seedB) < endPositionH(seedA) || beginPositionV(seedB) < endPositionV(seedA))
        return false;

    // Distance is maximal distance, this corresponds to going the
    // distacen in the smaller distance with matches/mismatches and
    // the rest with indels.
    TThreshold distance = _max(beginPositionH(seedB) - endPositionH(seedA), beginPositionV(seedB) - endPositionV(seedA));
    // Compare distance with threshold.
    return distance <= maxGapSize;
}


// Returns true iff b can be Chaos chained to a where a is the one to
// the upper left, b the one to the lower right.
//
// TODO(holtgrew): Replace bandwidth with diagonalDistance.
template <typename TSeedSpec, typename TSeedConfig, typename TDistanceThreshold, typename TBandwidthThreshold>
inline bool
_seedsCombineable(Seed<TSeedSpec, TSeedConfig> const & seedA,
                  Seed<TSeedSpec, TSeedConfig> const & seedB,
                  TDistanceThreshold const & maxGapSize,
                  TBandwidthThreshold const & bandwidth,
                  Chaos const &)
{
    // b has to be right of a for the two seeds to be chainable.
    if (beginPositionH(seedB) < endPositionH(seedA) || beginPositionV(seedB) < endPositionV(seedA))
        return false;

    // The diagonal distance has to be smaller than the bandwidth.
    // TODO(holtgrew): s/beginDiagonal/getBeginDiagonal/
    TBandwidthThreshold diagonalDistance = _abs(endDiagonal(seedB) - beginDiagonal(seedA));
    if (diagonalDistance > bandwidth)
        return false;

    // Distance is maximal distance, this corresponds to going the
    // distance in the smaller distance with matches/mismatches and
    // the rest with indels.
    TDistanceThreshold distance = _max(beginPositionH(seedB) - endPositionH(seedA), beginPositionV(seedB) - endPositionV(seedA));
    // Compare distance with threshold.
    return distance <= maxGapSize;
}


// Updating the coordinates of seeds is the same for merging and
// simple chaining.  Only the score computation differs.
template <typename TSeedConfig>
inline void
_updateSeedsCoordinatesMergeOrSimpleChain(
        Seed<Simple, TSeedConfig> & seed,
        Seed<Simple, TSeedConfig> const & other)
{
    setBeginPositionH(seed, std::min(beginPositionH(seed), beginPositionH(other)));
    setBeginPositionV(seed, std::min(beginPositionV(seed), beginPositionV(other)));
    setEndPositionH(seed, std::max(endPositionH(seed), endPositionH(other)));
    setEndPositionV(seed, std::max(endPositionV(seed), endPositionV(other)));
    setLowerDiagonal(seed, std::min(lowerDiagonal(seed), lowerDiagonal(other)));
    setUpperDiagonal(seed, std::max(upperDiagonal(seed), upperDiagonal(other)));
}


template <typename TSeedConfig, typename TScoreValue, typename TSequence0, typename TSequence1>
inline void
_combineSeeds(Seed<Simple, TSeedConfig> & seed,
              Seed<Simple, TSeedConfig> const & other,
              Score<TScoreValue, Simple> const & /*scoringScheme*/,
              TSequence0 const & /*sequence0*/,
              TSequence1 const & /*sequence1*/,
              Merge const &)
{
    _updateSeedsScoreMerge(seed, other);
    _updateSeedsCoordinatesMergeOrSimpleChain(seed, other);
}


template <typename TSeedConfig, typename TScoreValue, typename TSequence0, typename TSequence1>
inline void
_combineSeeds(Seed<Simple, TSeedConfig> & seed,
              Seed<Simple, TSeedConfig> const & other,
              Score<TScoreValue, Simple> const & scoringScheme,
              TSequence0 const & /*sequence0*/,
              TSequence1 const & /*sequence1*/,
              SimpleChain const &)
{
    _updateSeedsScoreSimpleChain(seed, other, scoringScheme);
    _updateSeedsCoordinatesMergeOrSimpleChain(seed, other);
}


template <typename TSeedConfig, typename TScoreValue, typename TSequence0, typename TSequence1>
inline void
_combineSeeds(Seed<Simple, TSeedConfig> & seed,
              Seed<Simple, TSeedConfig> const & other,
              Score<TScoreValue, Simple> const & scoringScheme,
              TSequence0 const & sequence0,
              TSequence1 const & sequence1,
              Chaos const &)
{
    typedef Seed<Simple, TSeedConfig> TSeed;
    typedef typename Position<TSeed>::Type TPosition;

    // TODO(holtgrew): Assert seed left of other.

    // Compute gaps in both dimensions, the remaining gap is the
    // vertical/horizontal distance we will not fill with CHAOS
    // chaining.
    //
    // TODO(holtgrew): We need + 1 here, do we need it anywhere else?
    TPosition gapDim0 = beginPositionH(other) - endPositionH(seed);
    TPosition gapDim1 = beginPositionV(other) - endPositionV(seed);
    TPosition minGap = _min(gapDim0, gapDim1);
    TPosition maxGap = _max(gapDim0, gapDim1);
    TPosition remainingGap = maxGap - minGap;

    // Compute new score using the CHAOS method.
    //
    // First, compute the score when force-aligning from seed.
    TPosition posLeft0 = endPositionH(seed);
    TPosition posLeft1 = endPositionV(seed);
    TScoreValue tmpScore = 0;
    // TODO(holtgrew): Probably better use iterators on sequences!
    for (TPosition i = 0; i < minGap; ++i)
        tmpScore += score(scoringScheme, sequence0[posLeft0 + i], sequence1[posLeft1 + i]);

    SEQAN_ASSERT_GT(beginPositionH(other), static_cast<TPosition>(0));
    SEQAN_ASSERT_GT(beginPositionV(other), static_cast<TPosition>(0));
    TPosition posRight0 = beginPositionH(other);
    TPosition posRight1 = beginPositionV(other);

    // Now, try to put the gap at each position and get the position
    // with the highest score.  If there are two such positions, the
    // first one found is returned which is the one that is furthest
    // away from seed.
    // TPosition bestGapPos = 0;  // delta to lowermost position  // TODO(holtgrew): Set but not used;  Remove?
    TScoreValue bestScore = tmpScore;
    for (TPosition i = 1; i < minGap; ++i) {
        tmpScore -= score(scoringScheme, sequence0[posLeft0 + minGap - i], sequence1[posLeft1 + minGap - i]);
        tmpScore += score(scoringScheme, sequence0[posRight0 - i], sequence1[posRight1 - i]);
        if (tmpScore > bestScore) {
            // Found a better score.
            bestScore = tmpScore;
            // bestGapPos = i;  // TODO(holtgrew): Set but not used;  Remove?
        }
    }

    // Now, the best gap is when extending the lower right seed
    // (other) by bestGapPos to the upper right.  However, this is
    // ignored for simple seeds: We simply update the score and are
    // done.
    _updateSeedsScoreChaos(seed, other, bestScore + remainingGap * scoreGap(scoringScheme));

    // For simple seeds, the coordinate computation is the same as for
    // merge/simple chain.
    //
    // TODO(holtgrew): Adjust the name of updateSeedsCoordinatesMergeOrSimpleChain to reflect this.
    _updateSeedsCoordinatesMergeOrSimpleChain(seed, other);
}


template <typename TSeedConfig, typename TScoreValue, typename TSequence0, typename TSequence1>
inline void
_combineSeeds(Seed<ChainedSeed, TSeedConfig> & seed,
              Seed<ChainedSeed, TSeedConfig> const & other,
              Score<TScoreValue, Simple> const & /*scoringScheme*/,
              TSequence0 const & /*sequence0*/,
              TSequence1 const & /*sequence1*/,
              Merge const &)
{
    // For chained seeds, we first remove all diagonals from seed
    // until the last diagonal of seed starts truly before other.
    // Then, we possibly shorten the last diagonal.  Finally, we copy
    // over all diagonals from other.

    // std::cout << "Merging chained seeds " << seed << " and " << other << std::endl;
    SEQAN_ASSERT_LEQ_MSG(beginPositionH(seed), beginPositionH(other), "Monotony in both dimensions required for merging.");
    SEQAN_ASSERT_LEQ_MSG(beginPositionV(seed), beginPositionV(other), "Monotony in both dimensions required for merging.");

    _updateSeedsScoreMerge(seed, other);

    // Remove diagonals.
    typedef Seed<ChainedSeed, TSeedConfig> TSeed;
    typedef typename Iterator<TSeed, Standard>::Type TIterator;
    TIterator it;
    // TODO(holtgrew): Could use back() instead of lastKept.
    TIterator lastKept = begin(seed);
    for (it = begin(seed); it != end(seed); ++it) {
        if (it->beginPositionH >= beginPositionH(other) && it->beginPositionV >= beginPositionV(other))
            break;
        lastKept = it;
    }
    if (it != end(seed))
        truncateDiagonals(seed, it);
    // std::cout << "Seed after truncating diagonals: " << seed << std::endl;

    // Shorten last diagonal if necessary.
    if (lastKept->beginPositionH + lastKept->length > beginPositionH(other) && lastKept->beginPositionV + lastKept->length > beginPositionV(other)) {
        lastKept->length = _min(beginPositionH(other) - lastKept->beginPositionH, beginPositionV(other) - lastKept->beginPositionV);
    } else if (lastKept->beginPositionH + lastKept->length > beginPositionH(other)) {
        lastKept->length = beginPositionH(other) - lastKept->beginPositionH;
    } else if (lastKept->beginPositionV + lastKept->length > beginPositionV(other)) {
        lastKept->length = beginPositionV(other) - lastKept->beginPositionV;
    }

    // Maybe remove shortened diagonal if its length is 0.
    if (back(seed).length == 0) {
        // TODO(holtgrew): Do not use dot method.
        seed._seedDiagonals.pop_back();
    }

    // Copy over other diagonals.
    typedef typename Iterator<TSeed const, Standard>::Type TConstIterator;
    for (TConstIterator it = begin(other, Standard()); it != end(other, Standard()); ++it)
        appendDiagonal(seed, *it);

    // std::cout << "Chained seed after merging: " << seed << std::endl;

    // TODO(holtgrew): Update lower and upper diagonals!
}


template <typename TSeedConfig, typename TScoreValue, typename TSequence0, typename TSequence1>
inline void
_combineSeeds(Seed<ChainedSeed, TSeedConfig> & seed,
              Seed<ChainedSeed, TSeedConfig> const & other,
              Score<TScoreValue, Simple> const & scoringScheme,
              TSequence0 const & /*sequence0*/,
              TSequence1 const & /*sequence1*/,
              SimpleChain const &)
{
    // Simply copy over the diagonals of the seed (other) into the
    // left one (seed) after updating the score.

    _updateSeedsScoreSimpleChain(seed, other, scoringScheme);

    // Copy over other diagonals.
    typedef Seed<ChainedSeed, TSeedConfig> TSeed;
    typedef typename Iterator<TSeed const, Standard>::Type TConstIterator;
    for (TConstIterator it = begin(other, Standard()); it != end(other, Standard()); ++it)
        appendDiagonal(seed, *it);
}


template <typename TSeedConfig, typename TScoreValue, typename TSequence0, typename TSequence1>
inline void
_combineSeeds(Seed<ChainedSeed, TSeedConfig> & seed,
              Seed<ChainedSeed, TSeedConfig> const & other,
              Score<TScoreValue, Simple> const & scoringScheme,
              TSequence0 const & sequence0,
              TSequence1 const & sequence1,
              Chaos const &)
{
    typedef Seed<ChainedSeed, TSeedConfig> TSeed;
    typedef typename Position<TSeed>::Type TPosition;
    typedef typename Iterator<TSeed const, Standard>::Type TConstIterator;

    // TODO(holtgrew): Assert seed left of other.

    // Compute gaps in both dimensions, the remaining gap is the
    // vertical/horizontal distance we will not fill with CHAOS
    // chaining.
    //
    // TODO(holtgrew): We need + 1 here, do we need it anywhere else?
    TPosition gapDim0 = beginPositionH(other) - endPositionH(seed);
    TPosition gapDim1 = beginPositionV(other) - endPositionV(seed);
    TPosition minGap = _min(gapDim0, gapDim1);
    TPosition maxGap = _max(gapDim0, gapDim1);
    TPosition remainingGap = maxGap - minGap;

    // Compute new score using the CHAOS method.
    //
    // First, compute the score when force-aligning from seed.
    TPosition posLeft0 = endPositionH(seed);
    TPosition posLeft1 = endPositionV(seed);
    TScoreValue tmpScore = 0;
    // TODO(holtgrew): Probably better use iterators on sequences!
    for (TPosition i = 0; i < minGap; ++i)
        tmpScore += score(scoringScheme, sequence0[posLeft0 + i], sequence1[posLeft1 + i]);

    SEQAN_ASSERT_GT(beginPositionH(other), static_cast<TPosition>(0));
    SEQAN_ASSERT_GT(beginPositionV(other), static_cast<TPosition>(0));
    TPosition posRight0 = beginPositionH(other);
    TPosition posRight1 = beginPositionV(other);

    // Now, try to put the gap at each position and get the position
    // with the highest score.  If there are two such positions, the
    // first one found is returned which is the one that is furthest
    // away from seed.
    TPosition bestGapPos = 0;  // delta to lowermost position
    TScoreValue bestScore = tmpScore;
    for (TPosition i = 1; i < minGap; ++i) {
        tmpScore -= score(scoringScheme, sequence0[posLeft0 + minGap - i], sequence1[posLeft1 + minGap - i]);
        tmpScore += score(scoringScheme, sequence0[posRight0 - i], sequence1[posRight1 - i]);
        if (tmpScore > bestScore) {
            // Found a better score.
            bestScore = tmpScore;
            bestGapPos = i;
        }
    }

    // Now, the best gap is when extending the lower right seed
    // (other) by bestGapPos to the upper right.  The upper left seed
    // is extended by (minGap - bestGapPos).
    //
    // Adjust last diagonal of seed.
    back(seed).length += minGap - bestGapPos;
    // Copy over the first diagonal of other and adjust diagonal.
    appendDiagonal(seed, front(other));
    back(seed).beginPositionH -= bestGapPos;
    back(seed).beginPositionV -= bestGapPos;
    back(seed).length += bestGapPos;
    // Copy over all other diagonals.
    TConstIterator it = begin(other, Standard());
    TConstIterator itEnd = end(other, Standard());
    // TODO(holtgrew): value(it) does not work here, the adaption around std::list needs more work!
    for (++it; it != itEnd; ++it)
        appendDiagonal(seed, *it);

    // Finally, we update the score and are done.
    _updateSeedsScoreChaos(seed, other, bestScore + remainingGap * scoreGap(scoringScheme));
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_SEEDS_COMBINATION_UNORDERED_H_
