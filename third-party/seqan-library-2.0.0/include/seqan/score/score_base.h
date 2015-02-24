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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================

// TODO(holtgrew): Should the public interface for the class Score not be defined here?

#ifndef SEQAN_SSCORE_BASE_H_
#define SEQAN_SSCORE_BASE_H_

namespace seqan {

/*!
 * @class Score
 * @headerfile <seqan/score.h>
 * @brief Scoring scheme.
 *
 * The Score class uses <b>similarity</b> scores, i.e. the greater the score value, the greater the similarity.
 * Depending on the exact score, the scores can also be negative for dissimilarity.  This choice blends in naturally
 * with the BLOSUM and PAM matrices, for examples.  For the edit distance common in computer science, this corresponds
 * to match scores of 0, mismatch and gap scores of -1.
 *
 * @signature template <[typename TValue[, typename TSpec]]>
 *            class Score;
 *
 * @tparam TValue The scoring value, defaults to <tt>int</tt>.
 * @tparam TSpec  The specialization type, defaults to <tt>Simple</tt>.
 *
 * @section Examples
 *
 * @include demos/score/score.cpp
 *
 * The output is as follows:
 *
 * @include demos/score/score.cpp.stdout
 */

template <typename TValue = int, typename TSpec = Simple>
class Score;

/*!
 * @mfn Score#Value
 * @brief Return the value type of the scoring scheme.
 *
 * @signature Value<TScore>::Type;
 *
 * @tparam TScore The Score specialization.
 *
 * @return Type The score value type of the scoring scheme.
 */

template <typename TValue, typename TSpec>
struct Value<Score<TValue, TSpec> > {
    typedef TValue Type;
};

// --------------------------------------------------------------------------
// Metafunction ConsensusScoreSequenceEntry
// --------------------------------------------------------------------------

/*!
 * @mfn Score#SequenceEntryForScore
 * @headerfile <seqan/score.h>
 *
 * @note This is used for unified interfaces for position dependent and independent scores.
 * @brief Returns representation type for a character of a position in a sequence.
 *
 * @signature SequenceEntryForScore<TScore, TSequence>::Type;
 *
 * @tparam TScore    The score type to use. Types: Score
 * @tparam TSequence The underlying sequence of the alignments or gaps. Types: ContainerConcept
 *
 * @return Type The type to use for representing a character in a sequence over its position.
 *
 *
 * @see Score#SequenceEntryForScore
 * @see Score#sequenceEntryForScore
 * @see ConsensusScoreSequenceEntry
 */

template <typename TScore, typename TSequence>
struct SequenceEntryForScore
{
    typedef typename Value<TSequence>::Type Type;
};

// --------------------------------------------------------------------------
// Function sequenceEntryForScore
// --------------------------------------------------------------------------

/*!
 * @fn Score#sequenceEntryForScore
 * @brief Helper function for element access, depending on score type.
 *
 * @signature TAlphabet sequenceEntryForScore(scoringScheme, seq, pos);
 *
 * @param[in] scoringScheme The Score to get the representation for.
 * @param[in] pos           The position of the character.
 * @param[in] seq           The sequence to get the representation for.
 *
 * @return TAlphabet The value of <tt>seq</tt> at <tt>pos</tt>.
 */

// TODO(rmaerker): Check if using iterator instead would be more efficient than subscript operator.
template <typename TScore, typename TSequence, typename TPosition>
inline typename Value<TSequence>::Type
sequenceEntryForScore(TScore const & /*scoringScheme*/, TSequence const & seq, TPosition pos)
{
    return seq[pos];
}

/*!
 * @fn Score#scoreGapOpenHorizontal
 * @brief Returns the score for opening a gap in horizontal direction.
 *
 * @signature TValue scoreGapOpenHorizontal(score, entryH, entryV);
 *
 * @param[in] score  The Score to query.
 * @param[in] entryH Entry in sequence one (horizontal), type from Score#SequenceEntryForScore.
 * @param[in] entryV Entry in sequence two (vertical), type from Score#SequenceEntryForScore.
 *
 * @return TValue The score open cost for gaps at the given position/entry.  TValue is the value type of score.
 *
 * @section Remarks
 *
 * Corresponds to a deletion event in sequence two and an insertion event in sequence one, respectively.
 */

template <typename TValue, typename TSpec, typename TSeqHValue, typename TSeqVValue>
inline TValue
scoreGapOpenHorizontal(
    Score<TValue, TSpec> const & me,
    TSeqHValue const & /*seqHVal*/,
    TSeqVValue const & /*seqHVal*/)
{
    SEQAN_CHECKPOINT;
    return scoreGapOpen(me);
}

/*!
 * @fn Score#scoreGapOpenVertical
 * @brief Returns the score for opening a gap in vertical direction.
 *
 * @signature TValue scoreGapOpenVertical(score, entryH, entryV);
 *
 * @param[in] score  The Score to query.
 * @param[in] entryH Entry in sequence one (horizontal), type from Score#SequenceEntryForScore.
 * @param[in] entryV Entry in sequence two (vertical), type from Score#SequenceEntryForScore.
 *
 * @return TValue The score open cost for gaps at the given position/entry.  TValue is the value type of score.
 *
 * @section Remarks
 *
 * Corresponds to a deletion event in sequence two and an insertion event in sequence one, respectively.
 */

template <typename TValue, typename TSpec, typename TSeqHValue, typename TSeqVValue>
inline TValue
scoreGapOpenVertical(
    Score<TValue, TSpec> const & me,
    TSeqHValue const & /*seqHVal*/,
    TSeqVValue const & /*seqHVal*/)
{
    SEQAN_CHECKPOINT;
    return scoreGapOpen(me);
}

/*!
 * @fn Score#scoreGapExtendHorizontal
 * @brief Returns the score for extending a gap in horizontal direction.
 *
 * @signature TValue scoreGapExtendHorizontal(score, entryH, entryV);
 *
 * @param[in] score  The Score to query.
 * @param[in] entryH Entry in sequence one (horizontal), type from Score#SequenceEntryForScore.
 * @param[in] entryV Entry in sequence two (vertical), type from Score#SequenceEntryForScore.
 *
 * @return TValue The score extension cost for gaps at the given position/entry.  TValue is the value type of score.
 *
 * @section Remarks
 *
 * Corresponds to a deletion event in sequence two and an insertion event in sequence one, respectively.
 */

template <typename TValue, typename TSpec, typename TSeqHValue, typename TSeqVValue>
inline TValue
scoreGapExtendHorizontal(
    Score<TValue, TSpec> const & me,
    TSeqHValue const & /*seqHVal*/,
    TSeqVValue const & /*seqHVal*/)
{
    SEQAN_CHECKPOINT;
    return scoreGapExtend(me);
}

/*!
 * @fn Score#scoreGapExtendVertical
 * @brief Returns the score for extending a gap in vertical direction.
 *
 * @signature TValue scoreGapExtendVertical(score, entryH, entryV);
 *
 * @param[in] score  The Score to query.
 * @param[in] entryH Entry in sequence one (horizontal), type from Score#SequenceEntryForScore.
 * @param[in] entryV Entry in sequence two (vertical), type from Score#SequenceEntryForScore.
 *
 * @return TValue The score extension cost for gaps at the given position/entry.  TValue is the value type of score.
 *
 * @section Remarks
 *
 * Corresponds to a deletion event in sequence one and an insertion event in sequence two, respectively.
 */

template <typename TValue, typename TSpec, typename TSeqHValue, typename TSeqVValue>
inline TValue
scoreGapExtendVertical(
    Score<TValue, TSpec> const & me,
    TSeqHValue const & /*seqHVal*/,
    TSeqVValue const & /*seqHVal*/)
{
    SEQAN_CHECKPOINT;
    return scoreGapExtend(me);
}

/*!
 * @fn Score#scoreGapHorizontal
 * @brief Returns the score for a gap in horizontal direction.
 *
 * @signature TValue scoreGapHorizontal(score, entryH, entryV);
 *
 * @param[in] score  The Score to query.
 * @param[in] entryH Entry in sequence one (horizontal), type from Score#SequenceEntryForScore.
 * @param[in] entryV Entry in sequence two (vertical), type from Score#SequenceEntryForScore.
 *
 * @return TValue The score gap cost for gaps at the given position/entry.  TValue is the value type of score.
 *
 * @section Remarks
 *
 * Corresponds to a deletion event in sequence two and an insertion event in sequence one, respectively.
 */

template <typename TValue, typename TSpec, typename TSeqHValue, typename TSeqVValue>
inline TValue
scoreGapHorizontal(
    Score<TValue, TSpec> const & me,
    TSeqHValue const & /*seqHVal*/,
    TSeqVValue const & /*seqHVal*/)
{
    SEQAN_CHECKPOINT;
    return scoreGap(me);
}

/*!
 * @fn Score#scoreGapVertical
 * @brief Returns the score for a gap in vertical direction.
 *
 * @signature TValue scoreGapVertical(score, entryH, entryV);
 *
 * @param[in] score  The Score to query.
 * @param[in] entryH Entry in sequence one (horizontal), type from Score#SequenceEntryForScore.
 * @param[in] entryV Entry in sequence two (vertical), type from Score#SequenceEntryForScore.
 *
 * @return TValue The score gap cost for gaps at the given position/entry.  TValue is the value type of score.
 *
 * @section Remarks
 *
 * Corresponds to a deletion event in sequence two and an insertion event in sequence one, respectively.
 */

template <typename TValue, typename TSpec, typename TSeqHValue, typename TSeqVValue>
inline TValue
scoreGapVertical(
    Score<TValue, TSpec> const & me,
    TSeqHValue const & /*seqHVal*/,
    TSeqVValue const & /*seqHVal*/)
{
    SEQAN_CHECKPOINT;
    return scoreGap(me);
}

/*!
 * @fn Score#score
 * @brief Returns similarity score for two sequence entries.
 *
 * @signature TValue score(score, entryH, entryV);
 *
 * @param[in] score  The Score to use for comparing the two sequence entries.
 * @param[in] entryH The entry in the first/horizontal sequence.
 * @param[in] entryV The entry in the second/vertical sequence.
 *
 * @return TValue The score similarity cost for gaps at the given position/entry.  TValue is the value type of
 *                score.
 */

template <typename TValue, typename TSpec, typename TSeqHVal, typename TSeqVVal>
inline TValue
score(Score<TValue, TSpec> const & me, TSeqHVal valH, TSeqVVal valV) {
    SEQAN_CHECKPOINT;
    if (valH == valV)
        return scoreMatch(me);
    else
        return scoreMismatch(me);
}

}  // namespace SEQAN_NAMESPACE_MAIN

#endif  // SEQAN_SSCORE_BASE_H_
