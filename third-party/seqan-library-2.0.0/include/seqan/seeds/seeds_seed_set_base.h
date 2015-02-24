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
// The class SeedSet.  ScoringScheme and related tags are defined in
// seeds_scoring_scheme.h
// ==========================================================================

#ifndef SEQAN_SEEDS_SEEDS_SEED_SET_BASE_H_
#define SEQAN_SEEDS_SEEDS_SEED_SET_BASE_H_

namespace seqan {

// ===========================================================================
// Forwards
// ===========================================================================

// TODO(holtgrew): Add DiagonalSorted as default specialization.

struct Unordered_;
typedef Tag<Unordered_> Unordered;

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// ---------------------------------------------------------------------------
// Class SeedSet
// ---------------------------------------------------------------------------

/*!
 * @class SeedSet
 * @headerfile <seqan/seeds.h>
 * @implements ContainerConcept
 * @brief Handles a set of seeds with local chaining on adding seeds.
 * @note At the moment only <tt>Unordered SeedSets</tt> are supported.
 *
 * @signature template <typename TSeed[, typename TSpec]>
 *            class SeedSet;
 *
 * @tparam TSeed     Type of the @link Seed @endlink objects stored in the seed set.
 * @tparam TSpec     Optional tag for seed set specialization. Defaults to <tt>Unordered</tt>.
 */

// ..param.TScored:Either UnScored or a seed set scoring scheme specification.
// ..param.TSpec:Specialization of the seed set.
// ..param.TSeedConfig:Configuration for the seeds.  Sensible defaults are chosen based on the other template parameters.

template <typename TSeed, typename TSpec = Unordered>
class SeedSet;

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

// Basic Container Functions

// TODO(holtgrew): dddoc {begin,end,length,front,back}All(T)

// SeedSet Functions

/*!
 * @fn SeedSet#addSeed
 *
 * @headerfile <seqan/seeds.h>
 *
 * @brief Adds a seed to an existing @link SeedSet @endlink using different
 *        algorithms for local chaining.
 *
 * @signature bool addSeed(seedSet, seed, distance, bandwidth, score, seqH, seqV, tag);
 * @signature bool addSeed(seedSet, seed, distance, score, SimpleChain);
 * @signature bool addSeed(seedSet, seed, distance, Merge);
 * @signature bool addSeed(seedSet, seed, Single);
 *
 * @param[in,out] seedSet   The SeedSet to add the seed to.
 * @param[in]     seed      The seed to be added.
 * @param[in]     distance  The maximal distance between the end point of the upper left and the begin point of the
 *                          lower right @link Seed @endlink allowed for local chaining.  NB: only Chaos, SimpleChain
 *                          and Merge require the distance information.
 * @param[in]     bandwidth The window size to search for a chainable @link Seed @endlink.  Note, only <tt>Chaos</tt>
 *                          requires the bandwidth information.
 * @param[in]     score     The scoring scheme.Note, only Chaos and SimpleChain require the score.
 *                          Type: @link SimpleScore @endlink.
 * @param[in]     seqH      Database sequence (horizontal).  Only required for Chaos Chaining.  Types: ContainerConcept.
 * @param[in]     seqV      Query sequence (vertical).  Only required for Chaos Chaining.  Types: ContainerConcept.
 * @param[in]     tag       Select the algorithm that is used to add the new seed.  Note that not every algorithm can
 *                          be used with each type of @link Seed @endlink.  See special signatures above.  The seed is
 *                          copied and then added.
 *
 * @return bool <tt>true</tt> if successful.  Adding can fail if no appropriate seed is there for chaining or merging.
 *              Adding using <tt>Single</tt> ever fails.
 *
 * @section Examples
 *
 * @include demos/seeds/seeds_add_seed.cpp
 *
 * The output is as follows:
 *
 * @code{.console}
 * Single Method:
 * Seed: Seed<Simple, TConfig>(4, 5, 8, 9, lower diag = -1, upper diag = -1)
 * Seed: Seed<Simple, TConfig>(10, 10, 15, 15, lower diag = 0, upper diag = 0)
 * Seed: Seed<Simple, TConfig>(14, 14, 18, 18, lower diag = 0, upper diag = 0)
 * Seed: Seed<Simple, TConfig>(21, 21, 24, 24, lower diag = 0, upper diag = 0)
 *
 *
 * Merge Method:
 * Seed: Seed<Simple, TConfig>(4, 5, 8, 9, lower diag = -1, upper diag = -1)
 * Seed: Seed<Simple, TConfig>(10, 10, 18, 18, lower diag = 0, upper diag = 0)
 * Seed: Seed<Simple, TConfig>(21, 21, 24, 24, lower diag = 0, upper diag = 0)
 *
 *
 * Chaos Method:
 * Seed: Seed<Simple, TConfig>(4, 5, 15, 15, lower diag = -1, upper diag = 0)
 * Seed: Seed<Simple, TConfig>(14, 14, 18, 18, lower diag = 0, upper diag = 0)
 * Seed: Seed<Simple, TConfig>(21, 21, 24, 24, lower diag = 0, upper diag = 0)
 * @endcode
 */


// ---------------------------------------------------------------------------
// Function minScore()
// ---------------------------------------------------------------------------

/*!
 * @fn SeedSet#minScore
 * @headerfile <seqan/seeds.h>
 * @brief Returns the threshold to distinguish between high-scoring and low-scoring seeds.
 *
 * @signature TSeedScore minScore(seedSet);
 *
 * @param[in] seedSet The SeedSet for which the threshold is set.  If the score of a seed is higher than the given
 *                    threshold, then it is virtually put into a container storing the high-scoring seeds which can
 *                    be iterated separately.
 *
 * @return TSeedScore The score threshold.  TSeedScore is the @link Seed#SeedScore @endlink of the contained seeds.
 *
 * @see SeedSet#setMinScore
 */

template <typename TSeed, typename TSeedSetSpec>
typename SeedScore<typename Value<SeedSet<TSeed, TSeedSetSpec> >::Type >::Type
minScore(SeedSet<TSeed, TSeedSetSpec> const & seedSet)
{
    return seedSet._minScore;
}

// ---------------------------------------------------------------------------
// Function setMinScore()
// ---------------------------------------------------------------------------

/*!
 * @fn SeedSet#setMinScore
 * @headerfile <seqan/seeds.h>
 * @brief Sets the threshold at which seeds are considered high-scoring.
 *
 * @signature void setMinScore(seedSet, scoreValue);
 *
 * @param[in,out] seedSet    The SeedSet for which the threshold is to be set.
 * @param[in]     scoreValue The new threshold to set.  If the score of a seed is higher than the given threshold, then
 *                           it is virtually put into a container storing the high-scoring seeds which can be iterated
 *                           separately  (@link IntegerConcept @endlink).
 *
 * @see SeedSet#minScore
 */

template <typename TSeed, typename TSeedSetSpec, typename TScoreValue>
void setMinScore(SeedSet<TSeed, TSeedSetSpec> & seedSet, TScoreValue val)
{
    seedSet._minScore = val;
}

// ---------------------------------------------------------------------------
// Function minSeedSize()
// ---------------------------------------------------------------------------

template <typename TSeed, typename TSeedSetSpec>
typename Size<typename Value<SeedSet<TSeed, TSeedSetSpec> >::Type >::Type
minSeedSize(SeedSet<TSeed, TSeedSetSpec> const & seedSet)
{
    return seedSet._minSeedSize;
}

// ---------------------------------------------------------------------------
// Function setMinSeedSize()
// ---------------------------------------------------------------------------

template <typename TSeed, typename TSeedSetSpec, typename TSize>
void setMinSeedSize(SeedSet<TSeed, TSeedSetSpec> & seedSet, TSize siz)
{
    seedSet._minSeedSize = siz;
}

// ---------------------------------------------------------------------------
// Helper Function _qualityReached()
// ---------------------------------------------------------------------------

// TODO(rmaerker): Is this function used anywhere?
template <typename TSeed, typename TSeedSetSpec, typename TSeedSpec, typename TSeedConfig>
inline bool _qualityReached(SeedSet<TSeed, TSeedSetSpec> const & seedSet,
                            Seed<TSeedSpec, TSeedConfig> const & seed)
{
    // TODO(rmaerker): If different seed configs are supported we must make sure, that the scoreValues are comparable to avoid compiler warnings.
    return score(seed) >= minScore(seedSet) && seedSize(seed) >= minSeedSize(seedSet);
}

// ---------------------------------------------------------------------------
// Function clear()
// ---------------------------------------------------------------------------

/*!
 * @fn SeedSet#clear
 * @headerfile <seqan/sequence.h>
 * @brief Clear the SeedSet.
 *
 * @signature void clear(seedSet);
 *
 * @param[in,out] seedSet The SeedSet to clear.
 */

template <typename TSeed, typename TSeedSetSpec>
inline void clear(SeedSet<TSeed, TSeedSetSpec> & seedSet)
{
    seedSet._seeds.clear();
    seedSet._minScore = 0;
    seedSet._minSeedSize = 0;
}


// Debugging / TikZ Output

template <typename TStream, typename TQuerySequence, typename TDatabaseSequence, typename TSeedSetSpec, typename TSeed>
inline void
__write(TStream & stream,
       TQuerySequence & sequence0,
       TDatabaseSequence & sequence1,
       SeedSet<TSeed, TSeedSetSpec> const & seedSet,
       Tikz_ const &)
{
//IOREV _nodoc_ specialization not documented
    typedef SeedSet<TSeed, TSeedSetSpec> TSeedSet;

    stream << "\\begin{tikzpicture}[" << std::endl
           << "    seed/.style={very thick}," << std::endl
           << "    seed diagonal/.style={red,<->}" << std::endl
           << "    ]" << std::endl;

    // Draw sequences.
    stream << "  \\draw";
    // Draw query / sequence 0;
    for (unsigned i = 0; i < length(sequence0); ++i)
        stream << std::endl << "    (0, -" << i << ") node {" << sequence0[i] << "}";
    stream << std::endl;
    // Draw database / sequence 1.
    for (unsigned i = 0; i < length(sequence1); ++i)
        stream << std::endl << "    (" << i << ", 0) node {" << sequence1[i] << "}";
    stream << ";" << std::endl;

    // Draw seeds.
    typedef typename Iterator<TSeedSet const, Standard>::Type TIterator;
    for (TIterator it = begin(seedSet); it != end(seedSet); ++it)
        _write(stream, value(it), Tikz_());
    stream << "\\end{tikzpicture}" << std::endl;
}

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_SEED_SET_BASE_H_
