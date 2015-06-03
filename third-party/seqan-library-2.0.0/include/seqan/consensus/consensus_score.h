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
// Author: Tobias Rausch <rausch@embl.de>
// ==========================================================================

#ifndef SEQAN_HEADER_SEQAN_CONSENSUS_SH
#define SEQAN_HEADER_SEQAN_CONSENSUS_SH


namespace SEQAN_NAMESPACE_MAIN
{


static const int SEQAN_CONSENSUS_UNITY = 1 << 20;

//////////////////////////////////////////////////////////////////////////////
// Consensus score tags
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

struct ConsensusScore_;
typedef Tag<ConsensusScore_> const ConsensusScore;

//////////////////////////////////////////////////////////////////////////////

struct FractionalScore_;
typedef Tag<FractionalScore_> const FractionalScore;

//////////////////////////////////////////////////////////////////////////////

template<typename TScore1, typename TScore2>
struct WeightedConsensusScore;




//////////////////////////////////////////////////////////////////////////////
// Scoring classes
//////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////
// ConsensusScore
//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
class Score<TValue, ConsensusScore>
{
public:
    String<TValue> consensus_set;        // Is the alphabet character part of the consensus set for each column

public:
    Score() {}

};

// --------------------------------------------------------------------------
// Class ConsensusScoreSequenceEntry
// --------------------------------------------------------------------------

/*!
 * @class ConsensusScoreSequenceEntry
 * @headerfile <seqan/consensus.h>
 * @brief Wrapper for a pointer to a sequence and a position in this sequence.
 *
 * @signature template <typename TSequence>
 *            class ConsensusScoreSequenceEntry;
 *
 * @tparam TSequence The sequence type this entry type is for.
 *
 * This is used for unified interfaces for position dependent and independetn scores.
 */

template <typename TSequence>
class ConsensusScoreSequenceEntry
{
public:
    typedef typename Position<TSequence>::Type TPosition;

    TSequence const * _seq;
    TPosition _pos;

    ConsensusScoreSequenceEntry() : _seq(), _pos(0)
    {}

    template <typename TPosition2>
    ConsensusScoreSequenceEntry(TSequence const & seq, TPosition2 pos) : _seq(&seq), _pos(pos)
    {}
};

/*!
 * @fn ConsensusScoreSequenceEntry#host
 * @brief Returns reference to sequence from entry.
 *
 * @signature TSequence host(entry);
 *
 * @param[in] entry The ConsensusScoreSequenceEntry to query for its host.
 *
 * @return TSequence A reference to the underlying sequence.
 */

template <typename TSequence>
inline TSequence const &
host(ConsensusScoreSequenceEntry<TSequence> & entry)
{
    return *entry._seq;
}

template <typename TSequence>
inline TSequence const &
host(ConsensusScoreSequenceEntry<TSequence> const & entry)
{
    return *entry._seq;
}

/*!
 * @fn ConsensusScoreSequenceEntry#position
 * @brief Returns position stored in <tt>entry</tt>.
 *
 * @signature TPosition position(entry);
 *
 * @param[in] entry The ConsensusScoreSequenceEntry to query.
 *
 * @return TPosition The position of the entry.  The type is
 *                   <tt>ConsensusScoreSequenceEntry&lt;TSequence&gt;::TPosition</tt>.
 */

template <typename TSequence>
inline typename ConsensusScoreSequenceEntry<TSequence>::TPosition
position(ConsensusScoreSequenceEntry<TSequence> const & entry)
{
    return entry._pos;
}

/*!
 * @fn ConsensusScoreSequenceEntry#value
 * @brief Returns value of character referenced by <tt>entry</tt>.
 *
 * @signature TValue value(entry);
 *
 * @param[in] entry The ConsensusScoreSequenceEntry to query.
 *
 * @return TValue The value of the sequence at the current position.
 */

template <typename TSequence>
inline typename Value<TSequence>::Type
value(ConsensusScoreSequenceEntry<TSequence> & entry)
{
    return host(entry)[position(entry)];
}

template <typename TSequence>
inline typename Value<TSequence>::Type
value(ConsensusScoreSequenceEntry<TSequence> const & entry)
{
    return host(entry)[position(entry)];
}

// --------------------------------------------------------------------------
// Metafunction SequenceEntryForScore                       [Consensus Score]
// --------------------------------------------------------------------------

template <typename TValue, typename TSequence>
struct SequenceEntryForScore<Score<TValue, ConsensusScore>, TSequence>
{
    typedef ConsensusScoreSequenceEntry<TSequence> Type;
};

template <typename TValue, typename TSequence>
struct SequenceEntryForScore<Score<TValue, ConsensusScore> const, TSequence> :
       SequenceEntryForScore<Score<TValue, ConsensusScore>, TSequence>
{};

// --------------------------------------------------------------------------
// Metafunction SequenceEntryForScore                      [Fractional Score]
// --------------------------------------------------------------------------

template <typename TValue, typename TSequence>
struct SequenceEntryForScore<Score<TValue, FractionalScore>, TSequence> :
       SequenceEntryForScore<Score<TValue, ConsensusScore>, TSequence>
{};

template <typename TValue, typename TSequence>
struct SequenceEntryForScore<Score<TValue, FractionalScore> const, TSequence> :
       SequenceEntryForScore<Score<TValue, ConsensusScore>, TSequence>
{};

// --------------------------------------------------------------------------
// Metafunction SequenceEntryForScore              [WeightedConsensus Score]
// --------------------------------------------------------------------------

template <typename TValue, typename TScore1, typename TScore2, typename TSequence>
struct SequenceEntryForScore<Score<TValue, WeightedConsensusScore<TScore1, TScore2> >, TSequence> :
       SequenceEntryForScore<Score<TValue, ConsensusScore>, TSequence>
{};

template <typename TValue, typename TScore1, typename TScore2, typename TSequence>
struct SequenceEntryForScore<Score<TValue, WeightedConsensusScore<TScore1, TScore2> > const, TSequence> :
       SequenceEntryForScore<Score<TValue, ConsensusScore>, TSequence>
{};

// --------------------------------------------------------------------------
// Function sequenceEntryForScore()                         [Consensus Score]
// --------------------------------------------------------------------------

template <typename TScoreValue, typename TSequence, typename TPosition>
inline ConsensusScoreSequenceEntry<TSequence>
sequenceEntryForScore(Score<TScoreValue, ConsensusScore> const &, TSequence const & seq, TPosition pos)
{
    return ConsensusScoreSequenceEntry<TSequence>(seq, pos);
}

template <typename TValue, typename TString>
inline void
assignProfile(Score<TValue, ConsensusScore>& me,
              TString const & profile)
{
//IOREV _notio_
    typedef typename Size<TString>::Type TSize;
    TSize alphSize = ValueSize<typename Value<TString>::Type>::VALUE;
    resize(me.consensus_set, alphSize * length(profile));

    typedef typename Iterator<TString const, Standard>::Type TIter;
    typedef typename Iterator<String<TValue>, Standard>::Type TConsSetIter;
    TConsSetIter itConsSet = begin(me.consensus_set, Standard());
    TIter it = begin(profile, Standard());
    TIter itEnd = end(profile, Standard());
    TSize maxCount = 0;
    for(;it!=itEnd;++it) {
        maxCount = 0;
        for(TSize i = 0; i<alphSize; ++i)
            if ((*it).count[i] > maxCount) maxCount = (*it).count[i];
        for(TSize i = 0; i<alphSize; ++i, ++itConsSet)
            *itConsSet = ((*it).count[i] == maxCount)? 0 : (-SEQAN_CONSENSUS_UNITY);
    }
}


template <typename TValue, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendHorizontal(
    Score<TValue, ConsensusScore> const & me,
    ConsensusScoreSequenceEntry<TSeq1> const & entry1,
    ConsensusScoreSequenceEntry<TSeq2> const & entry2)
{
    return ((int)position(entry2) < 0) ? -SEQAN_CONSENSUS_UNITY : me.consensus_set[position(entry1) * (ValueSize<typename Value<TSeq1>::Type>::VALUE) + (ValueSize<typename Value<TSeq1>::Type>::VALUE - 1)];
}

template <typename TValue, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenHorizontal(
    Score<TValue, ConsensusScore> const & me,
    ConsensusScoreSequenceEntry<TSeq1> const & entry1,
    ConsensusScoreSequenceEntry<TSeq2> const & entry2)
{
    return ((int)position(entry2) < 0) ? -2 * SEQAN_CONSENSUS_UNITY : 2 * me.consensus_set[position(entry1) * (ValueSize<typename Value<TSeq1>::Type>::VALUE) + (ValueSize<typename Value<TSeq1>::Type>::VALUE - 1)];
}

template <typename TValue, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendVertical(
    Score<TValue, ConsensusScore> const &,
    ConsensusScoreSequenceEntry<TSeq1> const & /*entry1*/,
    ConsensusScoreSequenceEntry<TSeq2> const & /*entry2*/)
{
    return -SEQAN_CONSENSUS_UNITY;
}

template <typename TValue, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenVertical(
    Score<TValue, ConsensusScore> const &,
    ConsensusScoreSequenceEntry<TSeq1> const & /*entry1*/,
    ConsensusScoreSequenceEntry<TSeq2> const & /*entry2*/)
{
    return -2 * SEQAN_CONSENSUS_UNITY;
}


template <typename TValue, typename TSeq1, typename TSeq2>
inline TValue
score(Score<TValue, ConsensusScore> const & me,
      ConsensusScoreSequenceEntry<TSeq1> const & entry1,
      ConsensusScoreSequenceEntry<TSeq2> const & entry2)
{
    return me.consensus_set[position(entry1) * (ValueSize<typename Value<TSeq1>::Type>::VALUE) + value(entry2).count[0]];
}









//////////////////////////////////////////////////////////////////////////////
// FractionalScore
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
class Score<TValue, FractionalScore>
{
public:
    String<int> sum;        // Total number of profile characters in each column

public:
    Score() {}
};


// --------------------------------------------------------------------------
// Function sequenceEntryForScore()                        [Fractional Score]
// --------------------------------------------------------------------------

template <typename TScoreValue, typename TSequence, typename TPosition>
inline ConsensusScoreSequenceEntry<TSequence>
sequenceEntryForScore(Score<TScoreValue, FractionalScore> const &, TSequence const & seq, TPosition pos)
{
    return ConsensusScoreSequenceEntry<TSequence>(seq, pos);
}

template <typename TValue, typename TString>
inline void
assignProfile(Score<TValue, FractionalScore> & me,
              TString const & profile)
{
//IOREV _notio_
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


template <typename TValue, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendHorizontal(
    Score<TValue, FractionalScore> const & me,
    ConsensusScoreSequenceEntry<TSeq1> const & entry1,
    ConsensusScoreSequenceEntry<TSeq2> const & entry2)
{
    return (( (int)position(entry2) < 0) || (!me.sum[position(entry1)])) ? -SEQAN_CONSENSUS_UNITY : ((TValue) (( (int)value(entry1).count[ValueSize<typename Value<TSeq1>::Type>::VALUE - 1] - me.sum[position(entry1)]) * SEQAN_CONSENSUS_UNITY) / me.sum[position(entry1)]);
}

template <typename TValue, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenHorizontal(
    Score<TValue, FractionalScore> const & me,
    ConsensusScoreSequenceEntry<TSeq1> const & entry1,
    ConsensusScoreSequenceEntry<TSeq2> const & entry2)
{
    return (( (int)position(entry2) < 0) || (!me.sum[position(entry1)])) ? -SEQAN_CONSENSUS_UNITY : ((TValue) (( (int) value(entry1).count[ValueSize<typename Value<TSeq1>::Type>::VALUE - 1] - me.sum[position(entry1)]) * SEQAN_CONSENSUS_UNITY) / me.sum[position(entry1)]);
}

template <typename TValue, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendVertical(
    Score<TValue, FractionalScore> const &,
    ConsensusScoreSequenceEntry<TSeq1> const & /*entry1*/,
    ConsensusScoreSequenceEntry<TSeq2> const & /*entry2*/)
{
    return -SEQAN_CONSENSUS_UNITY;
}

template <typename TValue, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenVertical(
    Score<TValue, FractionalScore> const &,
    ConsensusScoreSequenceEntry<TSeq1> const & /*entry1*/,
    ConsensusScoreSequenceEntry<TSeq2> const & /*entry2*/)
{
    return -SEQAN_CONSENSUS_UNITY;
}

template <typename TValue, typename TSeq1, typename TSeq2>
inline TValue
score(Score<TValue, FractionalScore> const & me,
      ConsensusScoreSequenceEntry<TSeq1> const & entry1,
      ConsensusScoreSequenceEntry<TSeq2> const & entry2)
{
    return (!me.sum[position(entry1)]) ? -SEQAN_CONSENSUS_UNITY : ((TValue) (((int)value(entry1).count[value(entry2).count[0]] - me.sum[position(entry1)]) * SEQAN_CONSENSUS_UNITY) / me.sum[position(entry1)]);
}








//////////////////////////////////////////////////////////////////////////////
// WeightedConsensusScore
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TScore1, typename TScore2>
class Score<TValue, WeightedConsensusScore<TScore1, TScore2> >
{
public:
    TScore1 sc1;
    TScore2 sc2;

public:
    Score() : sc1(), sc2() {}

};

// --------------------------------------------------------------------------
// Function sequenceEntryForScore()                 [WeightedConsensus Score]
// --------------------------------------------------------------------------

template <typename TScoreValue, typename TScore1, typename TScore2, typename TSequence, typename TPosition>
inline ConsensusScoreSequenceEntry<TSequence>
sequenceEntryForScore(Score<TScoreValue, WeightedConsensusScore<TScore1, TScore2> > const &,
                      TSequence const & seq,
                      TPosition pos)
{
    return ConsensusScoreSequenceEntry<TSequence>(seq, pos);
}

template <typename TValue, typename TScore1, typename TScore2, typename TString>
inline void
assignProfile(Score<TValue, WeightedConsensusScore<TScore1, TScore2> >& me,
              TString const& profile)
{
//IOREV _notio_
    assignProfile(me.sc1, profile);
    assignProfile(me.sc2, profile);
}



template <typename TValue, typename TScore1, typename TScore2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendHorizontal(
    Score<TValue, WeightedConsensusScore<TScore1, TScore2> > const & me,
    ConsensusScoreSequenceEntry<TSeq1> const & entry1,
    ConsensusScoreSequenceEntry<TSeq2> const & entry2)
{
    return (scoreGapExtendHorizontal(me.sc1, entry1, entry2) + scoreGapExtendHorizontal(me.sc2, entry1, entry2)) / (TValue) 2;
}

template <typename TValue, typename TScore1, typename TScore2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenHorizontal(
    Score<TValue, WeightedConsensusScore<TScore1, TScore2> > const & me,
    ConsensusScoreSequenceEntry<TSeq1> const & entry1,
    ConsensusScoreSequenceEntry<TSeq2> const & entry2)
{
    return (scoreGapOpenHorizontal(me.sc1, entry1, entry2) + scoreGapOpenHorizontal(me.sc2, entry1, entry2)) / (TValue) 2;
}


template <typename TValue, typename TScore1, typename TScore2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendVertical(
    Score<TValue, WeightedConsensusScore<TScore1, TScore2> > const & me,
    ConsensusScoreSequenceEntry<TSeq1> const & entry1,
    ConsensusScoreSequenceEntry<TSeq2> const & entry2)
{
    return (scoreGapExtendVertical(me.sc1, entry1, entry2) + scoreGapExtendVertical(me.sc2, entry1, entry2)) / (TValue) 2;
}

template <typename TValue, typename TScore1, typename TScore2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenVertical(
    Score<TValue, WeightedConsensusScore<TScore1, TScore2> > const & me,
    ConsensusScoreSequenceEntry<TSeq1> const & entry1,
    ConsensusScoreSequenceEntry<TSeq2> const & entry2)
{
    return (scoreGapOpenVertical(me.sc1, entry1, entry2) + scoreGapOpenVertical(me.sc2, entry1, entry2)) / (TValue) 2;
}


template <typename TValue, typename TScore1, typename TScore2, typename TSeq1, typename TSeq2>
inline TValue
score(Score<TValue, WeightedConsensusScore<TScore1, TScore2> > const & me,
      ConsensusScoreSequenceEntry<TSeq1> const & entry1,
      ConsensusScoreSequenceEntry<TSeq2> const & entry2)
{
    return (score(me.sc1, entry1, entry2) + score(me.sc2, entry1, entry2)) / (TValue) 2;
}

}

#endif

