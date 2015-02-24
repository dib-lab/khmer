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

#ifndef INCLUDE_SEQAN_CONSENSUS_OVERLAPPER_H_
#define INCLUDE_SEQAN_CONSENSUS_OVERLAPPER_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class OverlapperOptions_
// ----------------------------------------------------------------------------

// Overlapper configuration.
struct OverlapperOptions_
{
    OverlapperOptions_() : overlapErrorRate(0.05), overlapMinLength(40), logging(false)
    {}

    double overlapErrorRate;
    int overlapMinLength;
    bool logging;
};

// ----------------------------------------------------------------------------
// Class OverlapCandidate_
// ----------------------------------------------------------------------------

class OverlapCandidate_
{
public:
    unsigned seq0;
    unsigned seq1;
    int lDiag;
    int uDiag;

    // Constructors.
    OverlapCandidate_() : seq0(0), seq1(0), lDiag(0), uDiag(0) {}
    OverlapCandidate_(unsigned seq0, unsigned seq1, int lDiag, int uDiag) :
            seq0(seq0), seq1(seq1), lDiag(lDiag), uDiag(uDiag)
    {}

    // Compares two OverlapCanidate objects lexicographically.
    bool operator<(OverlapCandidate_ const & other) const
    { return std::make_pair(seq0, seq1) < std::make_pair(other.seq0, other.seq1); }
};

inline std::ostream & operator<<(std::ostream & out, OverlapCandidate_ const & cand)
{
    return out << "OverlapCandidate_(seq0=" << cand.seq0 << ", seq1=" << cand.seq1
               << ", lDiag=" << cand.lDiag << ", uDiag=" << cand.uDiag << ")";
}

// ----------------------------------------------------------------------------
// Class Overlap_
// ----------------------------------------------------------------------------

// Stores an overlap, usually such that seq0 is left of seq1 or if both start at the same position then seq0 is the id
// of the longer sequence (container).  In case of stacking, seq0 < seq1.

class Overlap_
{
public:
    // Value for invalid entry below.
    static const unsigned INVALID = (unsigned)-1;

    // Identifiers of the sequence.
    unsigned seq0;
    unsigned seq1;
    // Lengths of the sequences.
    unsigned len0;
    unsigned len1;
    // Begin positions of the sequences.
    unsigned begin0;
    unsigned begin1;
    // Edit distance errors of the alignment.
    unsigned errors;

    // Default constructor and construct with all values.
    Overlap_() : seq0(INVALID), seq1(INVALID), len0(INVALID), len1(INVALID), begin0(INVALID), begin1(INVALID),
                 errors(INVALID)
    {}

    Overlap_(unsigned seq0, unsigned seq1, unsigned len0, unsigned len1, unsigned begin0, unsigned begin1,
             unsigned errors = INVALID) :
            seq0(seq0), seq1(seq1), len0(len0), len1(len1), begin0(begin0), begin1(begin1), errors(errors)
    {}

    // Returns a "flipped" alignment, i.e. seq0 and seq1 change roles.
    Overlap_ flip() const { return Overlap_(seq1, seq0, len1, len0, begin1, begin0, errors); }

    // Returns a normalized overlap, i.e. if begin0 != begin1 then seq0 < seq1.  If begin0 == begin1 then ties are
    // broken by length (such that the container comes first), in case of stacks, ties are broken by id.
    Overlap_ normalize() const
    {
        if (begin0 != begin1)
        {
            return (begin0 < begin1) ? *this : flip();
        }
        else  // begin0 == begin1
        {
            if (len0 != len1)
                return (len0 > len1) ? *this : flip();
            else  // begin0 == begin1 && len0 == len1
                return (seq0 < seq1) ? *this : flip();
        }
    }

    // Returns the length of the overlap.
    int length() const
    {
        Overlap_ o = normalize();
        return (o.len0 - o.begin1);
    }

    // Returns "key", i.e. (seq0, seq1).
    std::pair<unsigned, unsigned> key() const { return std::make_pair(seq0, seq1); }

    // Compares two Overlap_ objects lexicographically.
    bool operator<(Overlap_ const & other) const { return (key() < other.key()); }
};

inline std::ostream & operator<<(std::ostream & out, Overlap_ const & ovl)
{
    return out << "Overlap_(seq0=" << ovl.seq0 << ", seq1=" << ovl.seq1 << ", len0=" << ovl.len0
               << ", len1=" << ovl.len1 << ", begin0=" << ovl.begin0 << ", begin1=" << ovl.begin1
               << ", errors=" << ovl.errors << ")";
}

// ----------------------------------------------------------------------------
// Class Overlapper_
// ----------------------------------------------------------------------------

template <typename TFragments, typename TSequence>
class Overlapper_
{
public:
    explicit Overlapper_(OverlapperOptions_ options = OverlapperOptions_()) : options(options)
    {}

    // Compute overlap from candidate and store it in overlap and alignment on success.  Return false if the overlap did
    // not exist.
    bool computeOverlap(Overlap_ & overlap,
                        TFragments & alignment,
                        TSequence const & seqH,
                        TSequence const & seqV,
                        OverlapCandidate_ const & candidate) const;

private:
    // Generate Overlap_ record from the alignment stored in fragments.  Length information is taken from seqs.
    Overlap_ overlapFromAlignment(
            String<Fragment<> > const & fragments,
            StringSet<TSequence, Dependent<> > const & strings) const;

    template <typename TSequenceH, typename TSequenceV, typename TAlignConfig, typename TAlgoTag>
    void _fixBandSize(int & lDiag,
                      int & uDiag,
                      TSequenceH const & seqH,
                      TSequenceV const & seqV,
                      TAlignConfig const & /*alignConfig*/,
                      TAlgoTag const & /*algoTag*/) const
    {
        // typedef typename SubstituteAlignConfig_<TAlignConfig>::Type TFreeEndGaps;
        typedef typename If<typename IsSameType<TAlgoTag, Gotoh>::Type, AffineGaps, LinearGaps>::Type TGapsType;
        typedef typename SetupAlignmentProfile_<DPGlobal, TAlignConfig, TGapsType, TracebackConfig_<SingleTrace, GapsLeft> >::Type TDPProfile;

        if (uDiag < -(int)length(seqV))
            uDiag = -(int)length(seqV);
        if (lDiag > (int)length(seqH))
            lDiag = length(seqV);

        if (uDiag < 0 && !IsFreeEndGap_<TDPProfile, DPFirstColumn>::VALUE)
            uDiag = 0;

        if (lDiag > 0 && !IsFreeEndGap_<TDPProfile, DPFirstRow>::VALUE)
            lDiag = 0;

        if (uDiag + (int)length(seqV) < (int)length(seqH) && !IsFreeEndGap_<TDPProfile, DPLastRow>::VALUE)
            uDiag = (int)length(seqH) - (int)length(seqV);

        if (lDiag + (int)length(seqV) > (int)length(seqH) && !IsFreeEndGap_<TDPProfile, DPLastColumn>::VALUE)
            lDiag = (int)length(seqH) - (int)length(seqV);
    }

    OverlapperOptions_ options;
};

template <typename TFragments, typename TSequence>
Overlap_ Overlapper_<TFragments, TSequence>::overlapFromAlignment(
        String<Fragment<> > const & fragments,
        StringSet<TSequence, Dependent<> > const & strings) const
{
    typedef StringSet<TSequence, Dependent<> > TStringSet;
    typedef typename Value<TStringSet>::Type TString;
    TStringSet & stringsNC = const_cast<StringSet<TSequence, Dependent<> > &>(strings);

    if (options.logging)
    {
        std::cerr << "FRAGMENTS\n";
        for (unsigned i = 0; i < length(fragments); ++i)
            std::cerr << "  Fragment(" << fragments[i].seqId1 << ", " << fragments[i].begin1
                      << ", " << fragments[i].seqId2 << ", " << fragments[i].begin2 << ", "
                      << fragments[i].len << ")\n";
    }

    // TODO(holtgrew): overlap length should actually be the length of the alignment

    TFragments frags = fragments;
    std::sort(begin(frags, seqan::Standard()), end(frags, seqan::Standard()));

    typename Value<TFragments>::Type frag0 = front(frags);  // first
    unsigned id0 = sequenceId(frag0, 0);
    unsigned id1 = sequenceId(frag0, 1);
    unsigned len0 = length(strings[0]);
    unsigned len1 = length(strings[1]);

    typedef int TPos;
    std::pair<TPos, TPos> range0(seqan::maxValue<TPos>(), seqan::minValue<TPos>());
    std::pair<TPos, TPos> range1 = range0;
    int errors = 0;
    typedef typename Iterator<TFragments, Standard>::Type TFragmentsIter;
    for (TFragmentsIter itF = begin(frags, seqan::Standard()); itF != end(frags, seqan::Standard()); ++itF)
    {
        // Get some shortcuts.
        TPos fLen = fragmentLength(*itF);
        TPos begin0 = fragmentBegin(*itF, id0);
        TPos begin1 = fragmentBegin(*itF, id1);

        // Count indels.
        if (itF != begin(frags, seqan::Standard()))
        {
            SEQAN_ASSERT_LEQ(range0.second, begin0);
            SEQAN_ASSERT_LEQ(range1.second, begin1);
            SEQAN_ASSERT_NEQ((begin0 != range0.second), (begin1 != range1.second));
            errors += (begin0 - range0.second);
            errors += (begin1 - range1.second);
        }

        // Update begin/end position in either read.
        range0.first = std::min(range0.first, begin0);
        range0.second = std::max(range0.second, begin0 + fLen);
        range1.first = std::min(range1.first, begin1);
        range1.second = std::max(range1.second, begin1 + fLen);

        // Count matches/mismatches.
        typedef typename Infix<TString>::Type TInfix;
        typedef typename Iterator<TInfix, Standard>::Type TInfixIter;
        TInfix label0 = label(*itF, stringsNC, id0);
        TInfix label1 = label(*itF, stringsNC, id1);
        for (TInfixIter it0 = begin(label0, seqan::Standard()), it1 = begin(label1, seqan::Standard());
             it0 != end(label0, seqan::Standard()); ++it0, ++it1)
            errors += ((seqan::Dna5)*it0 == 'N' ||
                       (seqan::Dna5)*it1 == 'N' ||
                       (seqan::Dna5)*it0 != (seqan::Dna5)*it1);
    }

    // In case that the alignment to the right aligns to a gap, flush left.
    int delta = std::min(range0.first, range1.first);
    range0.first -= delta;
    range1.first -= delta;

    SEQAN_ASSERT_MSG(range0.first == 0 || range1.first == 0, "One must start at beginning");
    // NB: Do not activate the following, does not have to be true, can end in alignment to gap.
    // SEQAN_ASSERT_MSG(range0.second == len0 || range1.second == len1, "One must end at last");

    TPos begin1 = range0.first, begin0 = range1.first;
    // int overlapLen = std::max(range0.second - range0.first, range1.second - range1.first);

    if (options.logging)
    {
        std::cerr << "range0 = (" << range0.first << ", " << range0.second << ")\n"
                  << "range1 = (" << range1.first << ", " << range1.second << ")\n";
    }

    return Overlap_(id0, id1, len0, len1, begin0, begin1, errors);
}

template <typename TFragments, typename TSequence>
inline bool Overlapper_<TFragments, TSequence>::computeOverlap(Overlap_ & overlap,
                                                               TFragments & frags,
                                                               TSequence const & seqH,
                                                               TSequence const & seqV,
                                                               OverlapCandidate_ const & candidate) const
{
    clear(frags);

    if (options.logging)
        std::cerr << "Computing overlap\n"
                  << "  seqH: " << seqH << "\n"
                  << "  seqV: " << seqV << "\n"
                  << "  cand: " << candidate << "\n";

    StringSet<TSequence, Dependent<> > pairSet;
    appendValue(pairSet, const_cast<TSequence &>(seqH));  // id 0
    appendValue(pairSet, const_cast<TSequence &>(seqV));  // id 1

    Score<int, Simple> scoringScheme(1000, -1000, -1001);

    AlignConfig<true, true, true, true> alignConfig;
    int uDiag = candidate.uDiag;
    int lDiag = candidate.lDiag;

    _fixBandSize(lDiag, uDiag, pairSet[0], pairSet[1], alignConfig, Gotoh());

    // DEBUG
    if (options.logging)
    {
        std::cerr << "\n\n(alignment of " << candidate << ")\n"
                  << "0:\t" << pairSet[0] << "\n"
                  << "1:\t" << pairSet[1] << "\n";
        typedef typename Value<StringSet<TSequence> >::Type TStringSeq;
        Align<TStringSeq> align;
        resize(rows(align), 2);
        setSource(row(align, 0), pairSet[0]);
        setSource(row(align, 1), pairSet[1]);
        globalAlignment(align, scoringScheme, alignConfig, lDiag, uDiag, Gotoh());
        std::cerr << "\n" << align << "\n";
    }
    // /DEBUG

    int overlapScore = globalAlignment(frags, pairSet, scoringScheme, alignConfig, lDiag, uDiag,
                                       NeedlemanWunsch());
    (void)overlapScore;
    if (empty(frags))
        return false;

    overlap = overlapFromAlignment(frags, pairSet);

    // Replace ids in fragments and overlap.
    overlap.seq0 = candidate.seq0;
    overlap.seq1 = candidate.seq1;
    for (unsigned i = 0; i < length(frags); ++i)
    {
        sequenceId(frags[i], 0) = candidate.seq0;
        sequenceId(frags[i], 1) = candidate.seq1;
    }

    int ovlLen = overlap.length();
    bool ok = ((options.overlapMinLength < 0 || ovlLen >= options.overlapMinLength) &&
               ((options.overlapErrorRate < 0 || ovlLen == 0)
                || (1.0 * overlap.errors / ovlLen) <= options.overlapErrorRate + 0.00001));
    if (options.logging)
        std::cerr << "Resulting overlap:\t" << overlap << " (passes quality? " << ok << ")\n";
    return ok;
}

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_CONSENSUS_OVERLAPPER_H_
