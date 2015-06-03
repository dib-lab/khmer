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

#ifndef SEQAN_INCLUDE_PROFILE_SPROFILE_SEQ_H_
#define SEQAN_INCLUDE_PROFILE_SPROFILE_SEQ_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ProfileSeqScore

struct ProfileSeqScore_;
typedef Tag<ProfileSeqScore_> ProfileSeqScore;

template <typename TValue>
class Score<TValue, ProfileSeqScore>;

template <typename TValue, typename TString>
inline void
assignProfile(Score<TValue, ProfileSeqScore> & me, TString const & profile);

// ProfileSeqFracScore

struct ProfileSeqFracScore_;
typedef Tag<ProfileSeqFracScore_> ProfileSeqFracScore;

template <typename TValue>
class Score<TValue, ProfileSeqFracScore>;

template <typename TValue, typename TString>
inline void
assignProfile(Score<TValue, ProfileSeqFracScore> & me, TString const & profile);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class ProfileSeq Score
// ----------------------------------------------------------------------------

struct ProfileSeqScore_;
typedef Tag<ProfileSeqScore_> ProfileSeqScore;

/*!
 * @class ProfileSeqScore ProfileSeq Score
 * @extends Score
 * @headerfile <seqan/align_profile.h>
 * @brief Score for sequence-to-profile alignments.
 *
 * Using this class, you can align sequences to profiles.  The profile is assumed to be in the horizontal direction
 * (first row), the sequence in the vertical direction (second row).
 *
 * Scoring works as follows.
 *
 * The integer <tt>SEQAN_CONSENSUS_UNITY</tt> and fractions thereof are used to express scores.  Gap opens in the
 * profile are scored proportional to the number of gaps in the profile with two times unity, gap extends with one times
 * unity at the position.
 *
 * Gap opens in the sequence are scored with two times unity, gap extends in the sequence with unity.  Alignments of
 * profile characters to sequence characters are scored with the fraction of profile characters that match the sequence
 * characters times unity.
 *
 * @signature template <typename TValue>
 *            class Score<TValue, ProfileSeqScore>;
 *
 * @tparam TValue The integer type to use for representing scores.
 *
 * @section Examples
 *
 * The following example uses the ProfileSeq Score to align a sequence against a profile.
 * Note that we print the gap state for each position since profiles cannot be printed to one stdout character.
 *
 * @include demos/align_profile/profile_seq_score.cpp
 *
 * The output is as follows:
 *
 * @code{.output}
 * score value = -2097152
 * gaps in profile/sequence
 * pos G   S
 * 0   0   0
 * 1   1   0
 * 2   0   0
 * 3   1   0
 * 4   0   0
 * 5   0   0
 * @endcode
 */

/*!
 * @fn ProfileSeqScore::Score
 * @brief Constructor
 *
 * @signature Score::Score();
 * @signature Score::Score(profile);
 *
 * @param[in] profile The profile to copy from (@link AllocString @endlink of @link ProfileChar @endlink objects).
 *
 * When providing <tt>profile</tt>, the function @link ProfileSeqScore#assignProfile @endlink is automatically used to
 * assign the profile to this class.
 */

template <typename TValue>
class Score<TValue, ProfileSeqScore>
{
  public:
    // A table of position x (ord value) giving the counts of the characters at the given positions.
    String<int> consensusSet;

    Score() {}

    // Construct given a profile string.
    template <typename TProfile>
    explicit
    Score(TProfile const & profile)
    {
        assignProfile(*this, profile);
    }
};

// ----------------------------------------------------------------------------
// Class ProfileSeqFrac Score
// ----------------------------------------------------------------------------

struct ProfileSeqFracScore_;
typedef Tag<ProfileSeqFracScore_> ProfileSeqFracScore;

/*!
 * @class ProfileSeqFracScore ProfileSeqFrac Score
 * @extends Score
 * @headerfile <seqan/align_profile.h>
 * @brief Score for sequence-to-profile alignments.
 *
 * Using this class, you can align sequences to profiles.  The profile is assumed to be in the horizontal direction
 * (first row), the sequence in the vertical direction (second row).
 *
 * Scoring works as follows.
 *
 * The integer <tt>SEQAN_CONSENSUS_UNITY</tt> and fractions thereof are used to express scores.  Gap opens in the
 * profile are scored proportional to the number of gaps in the profile two times unity, gap extends are scored
 * proportional to the number of gaps in the profile at the position.  Gap opens in the sequence are scored with two
 * times unity, gap extends with one times unity.
 *
 * @signature template <typename TValue>
 *            class Score<TValue, ProfileSeqFracScore>;
 *
 * @tparam TValue The integer type to use for representing scores.
 *
 * @section Examples
 *
 * The following example uses the ProfileSeqFrac Score to align a sequence against a profile.  Note that we print the
 * gap state for each position since profiles cannot be printed to one stdout character.
 *
 * @include demos/align_profile/profile_seq_frac_score.cpp
 *
 * The output is as follows:
 *
 * @code
 * score value = -2097152
 * gaps in profile/sequence
 * pos G   S
 * 0   0   0
 * 1   1   0
 * 2   0   0
 * 3   1   0
 * 4   0   0
 * 5   0   0
 * @endcode
 */

/*!
 * @fn ProfileSeqFracScore::Score
 * @brief Constructor
 *
 * @signature Score::Score();
 * @signature Score::Score(profile);
 *
 * @param[in] profile The profile to copy from (@link AllocString @endlink of @link ProfileChar @endlink objects).
 *
 * When providing <tt>profile</tt>, the function @link ProfileSeqFracScore#assignProfile @endlink is automatically used to
 * assign the profile to this class.
 */


template <typename TValue>
class Score<TValue, ProfileSeqFracScore>
{
  public:
    // Total number of profile characters in each column
    String<int> sum;

    Score() {}

    // Construct given a profile string.
    template <typename TProfile>
    explicit
    Score(TProfile const & profile)
    {
        assignProfile(*this, profile);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction SequenceEntryForScore                      [ProfileSeq Score]
// --------------------------------------------------------------------------

// Returns the type that holds a sequence entry.  This is used for abstracting away the access to sequence characters.

template <typename TValue, typename TSequence>
struct SequenceEntryForScore<Score<TValue, ProfileSeqScore>, TSequence>
{
    typedef ConsensusScoreSequenceEntry<TSequence> Type;
};

template <typename TValue, typename TSequence>
struct SequenceEntryForScore<Score<TValue, ProfileSeqScore> const, TSequence> :
            SequenceEntryForScore<Score<TValue, ProfileSeqScore>, TSequence>
{};

// --------------------------------------------------------------------------
// Metafunction SequenceEntryForScore                  [ProfileSeqFrac Score]
// --------------------------------------------------------------------------

// Returns the type that holds a sequence entry.  This is used for abstracting away the access to sequence characters.

template <typename TValue, typename TSequence>
struct SequenceEntryForScore<Score<TValue, ProfileSeqFracScore>, TSequence>
{
    typedef ConsensusScoreSequenceEntry<TSequence> Type;
};

template <typename TValue, typename TSequence>
struct SequenceEntryForScore<Score<TValue, ProfileSeqFracScore> const, TSequence> :
            SequenceEntryForScore<Score<TValue, ProfileSeqFracScore>, TSequence>
{};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function sequenceEntryForScore()                        [ProfileSeq Score]
// --------------------------------------------------------------------------

template <typename TScoreValue, typename TSequence, typename TPosition>
inline ConsensusScoreSequenceEntry<TSequence>
sequenceEntryForScore(Score<TScoreValue, ProfileSeqScore> const & /*sScheme*/,
                      TSequence const & seq, TPosition pos)
{
    return ConsensusScoreSequenceEntry<TSequence>(seq, pos);
}

// --------------------------------------------------------------------------
// Function assignProfile()                                [ProfileSeq Score]
// --------------------------------------------------------------------------

/*!
 * @fn ProfileSeqScore#assignProfile
 * @brief Assign profile to ProfileSeqScore.
 *
 * @signature void assignProfile(score, profile);
 *
 * @param[out] score   The ProfileSeqScore object to assign the profile for.
 * @param[in]  profile The profile to assign to the score.  @link AllocString @endlink of @link ProfileChar @endlink.
 */

template <typename TValue, typename TString>
inline void
assignProfile(Score<TValue, ProfileSeqScore> & me,
              TString const & profile)
{
    typedef typename Size<TString>::Type TSize;
    TSize alphSize = ValueSize<typename Value<TString>::Type>::VALUE;
    resize(me.consensusSet, alphSize * length(profile));

    typedef typename Iterator<TString const, Standard>::Type TIter;
    typedef typename Iterator<String<TValue>, Standard>::Type TConsSetIter;
    TConsSetIter itConsSet = begin(me.consensusSet, Standard());
    TIter it = begin(profile, Standard());
    TIter itEnd = end(profile, Standard());
    TSize maxCount = 0;
    for (;it!=itEnd;++it)
    {
        maxCount = 0;
        for (TSize i = 0; i<alphSize; ++i)
            if ((TSize)(*it).count[i] > maxCount)
                maxCount = (*it).count[i];
        for (TSize i = 0; i<alphSize; ++i, ++itConsSet)
            *itConsSet = ((TSize)(*it).count[i] == maxCount)? 0 : (-SEQAN_CONSENSUS_UNITY);
    }
}

// --------------------------------------------------------------------------
// Function scoreGapExtendHorizontal()                     [ProfileSeq Score]
// --------------------------------------------------------------------------

template <typename TValue, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendHorizontal(
        Score<TValue, ProfileSeqScore> const & me,
        ConsensusScoreSequenceEntry<TSeq1> const & entry1,
        ConsensusScoreSequenceEntry<TSeq2> const & entry2)
{
    typedef typename Value<TSeq1>::Type TValue1;
    if ((int)position(entry2) < 0)
        return -SEQAN_CONSENSUS_UNITY;
    else
        return me.consensusSet[position(entry1) * (ValueSize<TValue1>::VALUE) + (ValueSize<TValue1>::VALUE - 1)];
}

// --------------------------------------------------------------------------
// Function scoreGapOpenHorizontal()                       [ProfileSeq Score]
// --------------------------------------------------------------------------

template <typename TValue, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenHorizontal(
        Score<TValue, ProfileSeqScore> const & me,
        ConsensusScoreSequenceEntry<TSeq1> const & entry1,
        ConsensusScoreSequenceEntry<TSeq2> const & entry2)
{
    typedef typename Value<TSeq1>::Type TValue1;
    if ((int)position(entry2) < 0)
        return -2 * SEQAN_CONSENSUS_UNITY;
    else
        return 2 * me.consensusSet[position(entry1) * (ValueSize<TValue1>::VALUE) + (ValueSize<TValue1>::VALUE - 1)];
}

// --------------------------------------------------------------------------
// Function scoreGapExtendVertical()                       [ProfileSeq Score]
// --------------------------------------------------------------------------

template <typename TValue, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenVertical(
        Score<TValue, ProfileSeqScore> const &,
        ConsensusScoreSequenceEntry<TSeq1> const & /*entry1*/,
        ConsensusScoreSequenceEntry<TSeq2> const & /*entry2*/)
{
    return -2 * SEQAN_CONSENSUS_UNITY;
}

// --------------------------------------------------------------------------
// Function scoreGapOpenVertical()                         [ProfileSeq Score]
// --------------------------------------------------------------------------

template <typename TValue, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendVertical(
        Score<TValue, ProfileSeqScore> const &,
        ConsensusScoreSequenceEntry<TSeq1> const & /*entry1*/,
        ConsensusScoreSequenceEntry<TSeq2> const & /*entry2*/)
{
    return -SEQAN_CONSENSUS_UNITY;
}

// --------------------------------------------------------------------------
// Function score()                                        [ProfileSeq Score]
// --------------------------------------------------------------------------

template <typename TValue, typename TSeq1, typename TSeq2>
inline TValue
score(Score<TValue, ProfileSeqScore> const & me,
      ConsensusScoreSequenceEntry<TSeq1> const & entry1,
      ConsensusScoreSequenceEntry<TSeq2> const & entry2)
{
    typedef typename Value<TSeq1>::Type TValue1;
    return me.consensusSet[position(entry1) * (ValueSize<TValue1>::VALUE) + ordValue(value(entry2))];
}

// --------------------------------------------------------------------------
// Function sequenceEntryForScore()                    [ProfileSeqFrac Score]
// --------------------------------------------------------------------------

template <typename TScoreValue, typename TSequence, typename TPosition>
inline ConsensusScoreSequenceEntry<TSequence>
sequenceEntryForScore(Score<TScoreValue, ProfileSeqFracScore> const & /*sScheme*/,
                      TSequence const & seq, TPosition pos)
{
    return ConsensusScoreSequenceEntry<TSequence>(seq, pos);
}

// --------------------------------------------------------------------------
// Function assignProfile()                            [ProfileSeqFrac Score]
// --------------------------------------------------------------------------

/*!
 * @fn ProfileSeqFracScore#assignProfile
 * @brief Assign profile to ProfileSeqFrac Score.
 *
 * @signature void assignProfile(score, profile);
 *
 * @param[out] score   The ProfileSeqScore object to assign the profile for.
 * @param[in]  profile The profile to assign to the score.  @link AllocString @endlink of @link ProfileChar @endlink.
 */

template <typename TValue, typename TString>
inline void
assignProfile(Score<TValue, ProfileSeqFracScore> & me,
              TString const & profile)
{
    typedef typename Size<TString>::Type TSize;
    resize(me.sum, length(profile));
    typedef typename Iterator<TString const, Standard>::Type TIter;
    typedef typename Iterator<String<int>, Standard>::Type TSumIter;
    TSumIter itSum = begin(me.sum, Standard());
    TIter it = begin(profile, Standard());
    TIter itEnd = end(profile, Standard());
    for (; it!=itEnd; ++it, ++itSum)
    {
        *itSum = 0;
        for (TSize i = 0; i < (TSize) ValueSize<typename Value<TString>::Type>::VALUE; ++i)
            *itSum += (*it).count[i];
    }
}

// --------------------------------------------------------------------------
// Function scoreGapExtendHorizontal()                 [ProfileSeqFrac Score]
// --------------------------------------------------------------------------

template <typename TValue, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendHorizontal(
        Score<TValue, ProfileSeqFracScore> const & me,
        ConsensusScoreSequenceEntry<TSeq1> const & entry1,
        ConsensusScoreSequenceEntry<TSeq2> const & entry2)
{
    if (((int)position(entry2) < 0) || (!me.sum[position(entry1)]))
        return -SEQAN_CONSENSUS_UNITY;
    else
        return ((TValue) (( (int)value(entry1).count[ValueSize<typename Value<TSeq1>::Type>::VALUE - 1] - me.sum[position(entry1)]) * SEQAN_CONSENSUS_UNITY) / me.sum[position(entry1)]);
}

// --------------------------------------------------------------------------
// Function scoreGapOpenHorizontal()                   [ProfileSeqFrac Score]
// --------------------------------------------------------------------------

template <typename TValue, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenHorizontal(
        Score<TValue, ProfileSeqFracScore> const & me,
        ConsensusScoreSequenceEntry<TSeq1> const & entry1,
        ConsensusScoreSequenceEntry<TSeq2> const & entry2)
{
    if (((int)position(entry2) < 0) || (!me.sum[position(entry1)]))
        return -SEQAN_CONSENSUS_UNITY;
    else
        return ((TValue) (((int)value(entry1).count[ValueSize<typename Value<TSeq1>::Type>::VALUE - 1] - me.sum[position(entry1)]) * SEQAN_CONSENSUS_UNITY) / me.sum[position(entry1)]);
}

// --------------------------------------------------------------------------
// Function scoreGapExtendVertical()                   [ProfileSeqFrac Score]
// --------------------------------------------------------------------------

template <typename TValue, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenVertical(
        Score<TValue, ProfileSeqFracScore> const &,
        ConsensusScoreSequenceEntry<TSeq1> const & /*entry1*/,
        ConsensusScoreSequenceEntry<TSeq2> const & /*entry2*/)
{
    return -SEQAN_CONSENSUS_UNITY;
}

// --------------------------------------------------------------------------
// Function scoreGapOpenVertical()                     [ProfileSeqFrac Score]
// --------------------------------------------------------------------------

template <typename TValue, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendVertical(
        Score<TValue, ProfileSeqFracScore> const &,
        ConsensusScoreSequenceEntry<TSeq1> const & /*entry1*/,
        ConsensusScoreSequenceEntry<TSeq2> const & /*entry2*/)
{
    return -SEQAN_CONSENSUS_UNITY;
}

// --------------------------------------------------------------------------
// Function score()                                    [ProfileSeqFrac Score]
// --------------------------------------------------------------------------

template <typename TValue, typename TSeq1, typename TSeq2>
inline TValue
score(Score<TValue, ProfileSeqFracScore> const & me,
      ConsensusScoreSequenceEntry<TSeq1> const & entry1,
      ConsensusScoreSequenceEntry<TSeq2> const & entry2)
{
    if (!me.sum[position(entry1)])
        return -SEQAN_CONSENSUS_UNITY;
    else
        return ((TValue) (((int)value(entry1).count[ordValue(value(entry2))] - me.sum[position(entry1)]) * SEQAN_CONSENSUS_UNITY) / me.sum[position(entry1)]);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_PROFILE_SPROFILE_SEQ_H_
