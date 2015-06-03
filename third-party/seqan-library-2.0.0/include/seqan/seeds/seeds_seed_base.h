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
// The class Seed.
// ==========================================================================

#ifndef SEQAN_SEEDS_SEEDS_SEED_BASE_H_
#define SEQAN_SEEDS_SEEDS_SEED_BASE_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// ---------------------------------------------------------------------------
// Class DefaultSeedConfig
// ---------------------------------------------------------------------------

/*!
 * @class DefaultSeedConfig
 * @headerfile <seqan/seeds.h>
 * @brief Default configuration for seeds without score.
 *
 * @signature struct DefaultSeedConfig.
 *
 * The default definition is as follows.  You use this as a blueprint for your own TConfig struct for a Seed or SeedSet
 * class.
 *
 * @code{.cpp}
 * struct DefaultSeedConfig
 * {
 *     typedef size_t TPosition;
 *     typedef size_t TSize;
 *     typedef MakeSigned_<size_t>::Type TDiagonal;
 *     typedef int TScoreValue;
 * };
 * @endcode
 *
 * @see Seed
 *
 * @typedef DefaultSeedConfig::TPosition
 * @brief Position type to use in the seed.
 * @signature size_t DefaultSeedConfig::TPosition;
 *
 * @typedef DefaultSeedConfig::TSize
 * @brief Size type to use in the seed.
 * @signature size_t DefaultSeedConfig::TSize;
 *
 * @typedef DefaultSeedConfig::TDiagonal
 * @brief Type to use for diagonals (signed version of TSize).
 * @signature (...) DefaultSeedConfig::TDiagonal;
 *
 * @typedef DefaultSeedConfig::TScoreValue
 * @brief Type to use for storing score values.
 * @signature int DefaultSeedConfig::TScoreValue;
 */

// Default configuration for seeds without score.
struct DefaultSeedConfig
{
    typedef size_t TPosition;
    typedef size_t TSize;
    typedef MakeSigned_<size_t>::Type TDiagonal;
    typedef int TScoreValue;
};

// ---------------------------------------------------------------------------
// Class Seed
// ---------------------------------------------------------------------------

/*!
 * @class Seed
 * @headerfile <seqan/seeds.h>
 * @brief Stores the start and end positions in the horizonal and vertical dimension.
 *
 * @signature template <typename TSpec[, typename TConfig]>
 *            class Seed;
 *
 * @tparam TSpec   The seed specialization type.
 * @tparam TConfig The configuration object to use for this seed, defaults to @link DefaultSeedConfig @endlink.
 *
 * @section Examples
 *
 * The following example shows the usage of three seed extension algorithms using the tags <tt>MaxExtend</tt>,
 * <tt>UnGappedXDrop</tt>, and <tt>GappedXDrop</tt>.
 *
 * @include demos/seeds/seeds_extension.cpp
 *
 * The output is as follows:
 *
 * @include demos/seeds/seeds_extension.cpp.stdout
 *
 * Here is an example for global seed chaining:
 *
 * @include demos/seeds/seeds_chaining.cpp
 */

template <typename TSpec, typename TConfig = DefaultSeedConfig>
class Seed;

// ===========================================================================
// Metafunctions
// ===========================================================================

// ---------------------------------------------------------------------------
// Metafunction Position
// ---------------------------------------------------------------------------

/*!
 * @mfn Seed#Position
 * @brief The position type of a Seed.
 *
 * @signature Position<TSeed>::Type;
 *
 * @tparam TSeed The Seed type to query.
 *
 * @return Type The position type of <tt>TSeed</tt>.
 */

template <typename TSpec, typename TConfig>
struct Position<Seed<TSpec, TConfig> >
{
    typedef typename TConfig::TPosition Type;
};

template <typename TSpec, typename TConfig>
struct Position<Seed<TSpec, TConfig> const> : Position<Seed<TSpec, TConfig> >
{};

// ---------------------------------------------------------------------------
// Metafunction Size
// ---------------------------------------------------------------------------

/*!
 * @mfn Seed#Size
 * @brief The size type of a Seed.
 *
 * @signature Size<TSeed>::Type;
 *
 * @tparam TSeed The Seed to query for its size type.
 *
 * @return Type The size type of <tt>TSeed</tt>.
 */

template <typename TSpec, typename TConfig>
struct Size<Seed<TSpec, TConfig> >
{
    typedef typename TConfig::TSize Type;
};

template <typename TSpec, typename TConfig>
struct Size<Seed<TSpec, TConfig> const> : Size<Seed<TSpec, TConfig> >
{};

// ---------------------------------------------------------------------------
// Metafunction Diagonal
// ---------------------------------------------------------------------------

/*!
 * @mfn Seed#Diagonal
 * @brief The diagonal type of a Seed.
 *
 * @signature Diagonal<TSeed>::Type;
 *
 * @tparam TSeed The Seed to query for its diagonal type.
 *
 * @return Type The diagonal type of <tt>TSeed</tt>.
 */

template <typename T>
struct Diagonal;

template <typename TSpec, typename TConfig>
struct Diagonal<Seed<TSpec, TConfig> >
{
    typedef typename TConfig::TDiagonal Type;
};

template <typename TSpec, typename TConfig>
struct Diagonal<Seed<TSpec, TConfig> const>
        : Diagonal<Seed<TSpec, TConfig> > {};

// ---------------------------------------------------------------------------
// Metafunction SeedScore
// ---------------------------------------------------------------------------

/*!
 * @mfn Seed#SeedScore
 * @brief The seed score type of a Seed.
 *
 * @signature SeedScore<TSeed>::Type;
 *
 * @tparam TSeed The Seed to query for its seed score type.
 *
 * @return Type The score type of <tt>TSeed</tt>.
 */

template <typename T>
struct SeedScore;

template <typename TSpec, typename TConfig>
struct SeedScore<Seed<TSpec, TConfig> >
{
    typedef typename TConfig::TScoreValue Type;
};

template <typename TSpec, typename TConfig>
struct SeedScore<Seed<TSpec, TConfig> const> : SeedScore<Seed<TSpec, TConfig> >
{};

// ===========================================================================
// Functions
// ===========================================================================

// TODO(holtgrew): COULD introduce {get,set}{Begin,End}(dim, value), but probably only necessary to make consistent with multi dimensional chaining interface.

// ---------------------------------------------------------------------------
// Function beginPositionH()
// ---------------------------------------------------------------------------

/*!
 * @fn Seed#beginPositionH
 * @brief Returns the begin position of the seed in the database (horizontal direction).
 *
 * @signature TPosition beginPositionH(seed);
 *
 * @param[in] seed The seed to query.
 *
 * @return TPosition The horizontal begin position of type @link Seed#Position @endlink.
 */

// ---------------------------------------------------------------------------
// Function endPositionH()
// ---------------------------------------------------------------------------

/*!
 * @fn Seed#endPositionH
 * @brief Returns the end position of the seed in the database (horizontal direction).
 *
 * @signature TPosition endPositionH(seed);
 *
 * @param[in] seed The seed to query.
 *
 * @return TPosition The horizontal end position of type @link Seed#Position @endlink.
 */


// ---------------------------------------------------------------------------
// Function setBeginPositionH()
// ---------------------------------------------------------------------------

/*!
 * @fn Seed#setBeginPositionH
 * @brief Sets the begin position of the seed in the database (horizontal direction).
 *
 * @signature void setBeginPositionH(seed, pos);
 *
 * @param[in,out] seed The Seed to set the horizontal begin position for.
 * @param[in]     pos  The value to set for the horizontal begin position.
 */

// ---------------------------------------------------------------------------
// Function setEndPositionH()
// ---------------------------------------------------------------------------

/*!
 * @fn Seed#setEndPositionH
 * @brief Sets the end position of the seed in the database (horizontal direction).
 *
 * @signature void setEndPositionH(seed, pos);
 *
 * @param[in,out] seed The Seed to set the horizontal end position for.
 * @param[in]     pos  The value to set for the horizontal end position.
 */

// ---------------------------------------------------------------------------
// Function beginPositionV()
// ---------------------------------------------------------------------------

/*!
 * @fn Seed#beginPositionV
 * @brief Returns the begin position of the seed in the query (vertical direction).
 *
 * @signature TPosition beginPositionV(seed);
 *
 * @param[in] seed The seed to query.
 *
 * @return TPosition The vertical begin position of type @link Seed#Position @endlink.
 */

// ---------------------------------------------------------------------------
// Function endPositionV()
// ---------------------------------------------------------------------------

/*!
 * @fn Seed#endPositionV
 * @brief Returns the end position of the seed in the query (vertical direction).
 *
 * @signature TPosition endPositionV(seed);
 *
 * @param[in] seed The seed to query.
 *
 * @return TPosition The vertical end position of type @link Seed#Position @endlink.
 */

// ---------------------------------------------------------------------------
// Function setBeginPositionV()
// ---------------------------------------------------------------------------

/*!
 * @fn Seed#setBeginPositionV
 * @brief Sets the begin position of the seed in the query (vertical direction).
 *
 * @signature void setBeginPositionV(seed, pos);
 *
 * @param[in,out] seed The Seed to set the vertical begin position for.
 * @param[in]     pos  The value to set for the vertical begin position.
 */

// ---------------------------------------------------------------------------
// Function setEndPositionV()
// ---------------------------------------------------------------------------

/*!
 * @fn Seed#setEndPositionV
 * @brief Sets the end position of the seed in the query (vertical direction).
 *
 * @signature void setEndPositionV(seed, pos);
 *
 * @param[in,out] seed The Seed to set the vertical end position for.
 * @param[in]     pos  The value to set for the vertical end position.
 */

// ---------------------------------------------------------------------------
// Function lowerDiagonal()
// ---------------------------------------------------------------------------

/*!
 * @fn Seed#lowerDiagonal
 * @brief Returns the leftmost diagonal of the seed (minimum diagonal value).
 *
 * @signature TDiagonal lowerDiagonal(seed);
 *
 * @param[in] seed The seed to query for its lower diagonal.
 *
 * @return TDiagonal The lower diagonal value of type @link Seed#Diagonal @endlink.
 */

template <typename TSpec, typename TConfig>
inline typename Diagonal<Seed<TSpec, TConfig> >::Type
lowerDiagonal(Seed<TSpec, TConfig> const & seed)
{
    return seed._lowerDiagonal;
}

// ---------------------------------------------------------------------------
// Function setLowerDiagonal()
// ---------------------------------------------------------------------------

/*!
 * @fn Seed#setLowerDiagonal
 * @brief Sets a new value for the leftmost diagonal of a seed.
 *
 * @signature void setLowerDigonal(seed, diag);
 *
 * @param[in,out] seed The Seed to set the diagonal value for.
 * @param[in]     diag The value to set for the diagonal.
 */

template <typename TSpec, typename TConfig, typename TDiagonal>
inline void
setLowerDiagonal(Seed<TSpec, TConfig> & seed, TDiagonal newDiag)
{
    seed._lowerDiagonal = newDiag;
}

// ---------------------------------------------------------------------------
// Function upperDiagonal()
// ---------------------------------------------------------------------------

/*!
 * @fn Seed#upperDiagonal
 * @brief Returns the rightmost diagonal of the seed (maximum diagonal value).
 *
 * @signature TDiagonal upperDiagonal(seed);
 *
 * @param[in] seed The seed to query for its upper diagonal.
 *
 * @return TDiagonal The upper diagonal value of type @link Seed#Diagonal @endlink.
 */

template <typename TSpec, typename TConfig>
inline typename Diagonal<Seed<TSpec, TConfig> >::Type
upperDiagonal(Seed<TSpec, TConfig> const & seed)
{
    return seed._upperDiagonal;
}

// ---------------------------------------------------------------------------
// Function setUpperDiagonal()
// ---------------------------------------------------------------------------

/*!
 * @fn Seed#setUpperDiagonal
 * @brief Sets a new value for the rightmost diagonal of a seed.
 *
 * @signature void setUpperDigonal(seed, diag);
 *
 * @param[in,out] seed The Seed to set the diagonal value for.
 * @param[in]     diag The value to set for the diagonal.
 */

template <typename TSpec, typename TConfig, typename TPosition>
inline void
setUpperDiagonal(Seed<TSpec, TConfig> & seed,
                 TPosition newDiag)
{
    seed._upperDiagonal = newDiag;
}

// ---------------------------------------------------------------------------
// Function seedSize()
// ---------------------------------------------------------------------------

/*!
 * @fn Seed#seedSize
 * @brief Returns the number of matches and mismatches of the seed.  This is the longest true diagonal fitting into
 *        its dimensions.
 *
 * @signature TSize seedSize(seed);
 *
 * @param[in] seed The Seed to query for its size.
 *
 * @return TSize The size of the type @link Seed#Size @endlink.
 */

template <typename TSpec, typename TConfig>
inline typename Size<Seed<TSpec, TConfig> >::Type
seedSize(Seed<TSpec, TConfig> & seed)
{
    return _max(endPositionH(seed) - beginPositionH(seed), endPositionV(seed) - beginPositionV(seed));
}

template <typename TSpec, typename TConfig>
inline typename Size<Seed<TSpec, TConfig> >::Type
seedSize(Seed<TSpec, TConfig> const & seed)
{
    return _max(endPositionH(seed) - beginPositionH(seed), endPositionV(seed) - beginPositionV(seed));
}

// ---------------------------------------------------------------------------
// Function beginDiagonal()
// ---------------------------------------------------------------------------

// Computed values, based on properties returned by getters.

/*!
 * @fn Seed#beginDiagonal
 * @brief Returns the begin diagonal of a Seed.
 *
 * @signature TDiagonal beginDiagonal(seed);
 *
 * @param[in] seed The Seed to query for its begin diagonal.
 *
 * @return TDiagonal The diagonal of the Seed's begin position of type @link Seed#Diagonal @endlink.
 */

// TODO(holtgrew): Rename to getBeginDiagonal.
template <typename TSpec, typename TConfig>
inline typename Diagonal<Seed<TSpec, TConfig> >::Type
beginDiagonal(Seed<TSpec, TConfig> const & seed)
{
    return beginPositionH(seed) - beginPositionV(seed);
}

// ---------------------------------------------------------------------------
// Function endDiagonal()
// ---------------------------------------------------------------------------

/*!
 * @fn Seed#endDiagonal
 * @brief Returns the end diagonal of a Seed.
 *
 * @signature TDiagonal endDiagonal(seed);
 *
 * @param[in] seed The Seed to query for its end diagonal.
 *
 * @return TDiagonal The diagonal of the Seed's end position of type @link Seed#Diagonal @endlink.
 */

template <typename TSpec, typename TConfig>
inline typename Diagonal<Seed<TSpec, TConfig> >::Type
endDiagonal(Seed<TSpec, TConfig> const & seed)
{
    return endPositionH(seed) - endPositionV(seed);
}

// ---------------------------------------------------------------------------
// Function score()
// ---------------------------------------------------------------------------

/*!
 * @fn Seed#score
 * @brief Returns the score of a Seed.
 *
 * @signature TSeedScore score(seed);
 *
 * @param[in] seed The Seed to query for its score.
 *
 * @return TSeedScore The score value of the seed of type @link Seed#SeedScore @endlink.
 */

template <typename TSeed>
inline typename SeedScore<TSeed>::Type
score(TSeed const & seed)
{
    return seed._score;
}

// ---------------------------------------------------------------------------
// Function setScore()
// ---------------------------------------------------------------------------

/*!
 * @fn Seed#setScore
 * @brief Sets the Seed score value.
 *
 * @signature void setScore(seed, scoreVal);
 *
 * @param[in,out] seed     The Seed to set the score value of.
 * @param[in]     scoreVal The score value to set.  The type can queried from the type of <tt>seed</tt> using
 *                         @link Seed#SeedScore @endlink.
 */

template <typename TSeed, typename TScore>
inline void
setScore(TSeed & seed, TScore const & score)
{
    seed._score = score;
}

// ---------------------------------------------------------------------------
// Functor LessBeginDiagonal_()
// ---------------------------------------------------------------------------

// Compare two seeds based on begin diagonal.

template <typename TSeed>
struct LessBeginDiagonal
{
    inline bool operator()(TSeed const & lhs, TSeed const & rhs) const
    {
        return beginDiagonal(lhs) < beginDiagonal(rhs);
    }
};

// ---------------------------------------------------------------------------
// Score Helper Functions
// ---------------------------------------------------------------------------

// Seed has score.  The assumption is that the score is proportional to the size of the seed.  Each seed contributes a
// fraction of its score that is proportional to the fraction of the score it contributes.
template <typename TSpec, typename TConfig>
inline void
_updateSeedsScoreMerge(Seed<TSpec, TConfig> & seed, Seed<TSpec, TConfig> const & other)
{
    typedef Seed<TSpec, TConfig> TSeed;
    typedef typename Size<TSeed>::Type TSize;

    // Compute new size.
    TSize newBegin0 = _min(beginPositionH(seed), beginPositionH(other));
    TSize newEnd0 = _max(endPositionH(seed), endPositionH(other));
    TSize newBegin1 = _min(beginPositionV(seed), beginPositionV(other));
    TSize newEnd1 = _max(endPositionV(seed), endPositionV(other));
    TSize newSize = _max(newEnd0 - newBegin0, newEnd1 - newBegin1);
    // New seed should be larger than either old one and overlap should be > 0.
    SEQAN_ASSERT_GEQ(newSize, seedSize(seed));
    SEQAN_ASSERT_GEQ(newSize, seedSize(other));
    SEQAN_ASSERT_LEQ(newSize, seedSize(seed) + seedSize(other));
    TSize overlap = seedSize(seed) + seedSize(other) - newSize;
    // Overlap should be smaller than or equal to either seed size.
    SEQAN_ASSERT_GEQ(seedSize(seed), overlap);
    SEQAN_ASSERT_GEQ(seedSize(other), overlap);

    // Compute fraction each seed contributes.
    TSize total = seedSize(seed) + seedSize(other) - overlap;
    double fracSeed = static_cast<double>(seedSize(seed) - 0.5 * overlap) / static_cast<double>(total);
    double fracOther = static_cast<double>(seedSize(other) - 0.5 * overlap) / static_cast<double>(total);
    typedef typename SeedScore<TSeed>::Type TScoreValue;
    TScoreValue newScore = static_cast<TScoreValue>(round(fracSeed * score(seed) + fracOther * score(other)));
    setScore(seed, newScore);
}

// Update the score of seed according to the gap size.
template <typename TSpec, typename TConfig, typename TScoreValue>
inline void
_updateSeedsScoreSimpleChain(Seed<TSpec, TConfig> & seed,
                             Seed<TSpec, TConfig> const & other,
                             Score<TScoreValue, Simple> const & scoringScheme)
{
    typedef Seed<TSpec, TConfig> TSeed;
    typedef typename Size<TSeed>::Type TSize;

    // TODO(holtgrew): seed must be the one to the upper left.

    // Only linear gap costs are supported.
    SEQAN_ASSERT_EQ(scoreGapOpen(scoringScheme), scoreGapExtend(scoringScheme));
    // Scores for gaps and mismatches must be penalties.
    SEQAN_ASSERT_LT(scoreGap(scoringScheme), static_cast<TScoreValue>(0));
    SEQAN_ASSERT_LT(scoreMismatch(scoringScheme), static_cast<TScoreValue>(0));

    // We use a simple heuristic for updating the seed's score the
    // systematically overestimates the gap costs: We close the gap
    // with a maximal diagonal and then remaining indels or just
    // indels, whatever yields a lower score.
    TSize maxDist = _max(beginPositionH(other) - endPositionH(seed), beginPositionV(other) - endPositionV(seed));
    TSize minDist = _max(beginPositionH(other) - endPositionH(seed), beginPositionV(other) - endPositionV(seed));
    TSize diagLen = minDist;
    TSize indelLen = maxDist - minDist;
    TScoreValue gapScore1 = diagLen * scoreMismatch(scoringScheme) + indelLen * scoreGap(scoringScheme);
    TScoreValue gapScore2 = (maxDist + minDist) * scoreGap(scoringScheme);
    TScoreValue gapScore = _max(gapScore1, gapScore2);

    // The new score is the sum of the seed scores and the better of
    // the gap scores computed above.
    setScore(seed, score(seed) + score(other) + gapScore);
}

// The new score is the sum of the other seed's score and the delta as computed by the chaining routine.
template <typename TSpec, typename TConfig, typename TScoreValue>
inline void
_updateSeedsScoreChaos(Seed<TSpec, TConfig> & seed, Seed<TSpec, TConfig> const & other, TScoreValue const & scoreDelta)
{
    setScore(seed, score(seed) + score(other) + scoreDelta);
}

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_SEED_BASE_H_
