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

// TODO(holtgrew): After this has been thoroughly tested, debug output should be removed.
// TODO(holtgrew): With gaps instead of align we can now create "real" profile-sequence alignments.
// TODO(holtgrew): For windows, we only need to create the profile of the window, this should save time for string resizing on replace().
// TODO(holtgrew): Need to resort after display, changes order and might break realignment!

#include <seqan/align_profile/score_profile_seq.h>

#ifndef INCLUDE_SEQAN_REALIGN_REALIGN_BASE_H_
#define INCLUDE_SEQAN_REALIGN_REALIGN_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class RealignmentOptions_
// ----------------------------------------------------------------------------

// Class used for holding the realignment options.

struct RealignmentOptions_
{
    // Enum for selecting a realignment method.
    enum RealignmentMethod
    {
        ANSON_MYERS_NW,    // Anson-Myers with Needleman-Wunsch alignment.
        ANSON_MYERS_GOTOH  // Anson-Mysers with Gotoh alignment.
    };

    // The realignment method to use.
    RealignmentMethod method;
    // The number of bases in the environment to use left and right of the current alignment.
    unsigned environment;
    // The bandwidth to use for realignment for lower and higher band.
    unsigned bandwidth;
    // Whether or not to include the reference as a read on the contig.  This is useful for disconnected read alignment
    // groups and gap-rich alignments (454 data).
    bool includeReference;
    // Whether or not to print timing information.
    bool printTiming;
    // Whether or not to print debugging information.
    bool debug;

    RealignmentOptions_() :
            method(ANSON_MYERS_NW), environment(0), bandwidth(0), includeReference(false), printTiming(false),
            debug(false)
    {}
};

// ----------------------------------------------------------------------------
// Class ContigAlignmentsInfo_
// ----------------------------------------------------------------------------

// Helper class that bundles the information for the alignments on a contig.

template <typename TFragmentStore>
struct ContigAlignmentsInfo_
{
    typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type TAlignedReadIter;
    typedef typename TFragmentStore::TReadPos TReadPos;

    TReadPos minPos, maxPos;
    TAlignedReadIter alignedItBegin, alignedItEnd;

    ContigAlignmentsInfo_() : minPos(0), maxPos(0), alignedItBegin(), alignedItEnd()
    {
        reset();
    }

    // Reset everything.
    void reset()
    {
        minPos = maxValue<TReadPos>();
        maxPos = minValue<TReadPos>();
        alignedItBegin = TAlignedReadIter();
        alignedItEnd = TAlignedReadIter();
    }
};

// ----------------------------------------------------------------------------
// Helper class for storing timing information.
// ----------------------------------------------------------------------------

struct AnsonMyersTimes_
{
    // Time for preparation step in realigner run() function.
    double mainBeginContig;
    // Time for writing back information in realigner run() function.
    double mainEndContig;
    // Total time spent in main.
    double mainTotal;

    // Time for extracting information from the store into the profile.
    double realignExtractProfile;
    // Time for realignment.
    double realignAlign;
    // Time for integration of information.
    double realignIntegration;
    // Time for total realignment rounds.
    double realignTotal;

    // The times for each realignment round.
    std::vector<double> roundTimes;

    AnsonMyersTimes_() :
            mainBeginContig(0), mainEndContig(0), mainTotal(0), realignExtractProfile(0),
            realignAlign(0), realignIntegration(0), realignTotal(0)
    {}
};


// ----------------------------------------------------------------------------
// Class ExtractProfileInfo_
// ----------------------------------------------------------------------------

// Struct for storing positions to extract from profile.

struct ExtractProfileInfo_
{
    // Begin and end positions in the profile to extract.
    unsigned profBeginPos, profEndPos;
    // Begin and end positions of the alignment within the profile.
    unsigned aliBeginPos, aliEndPos;

    ExtractProfileInfo_() : profBeginPos(0), profEndPos(0), aliBeginPos(0), aliEndPos(0)
    {}
};


// ----------------------------------------------------------------------------
// Class WindowClippingInfo_
// ----------------------------------------------------------------------------

// Information that stores flags for whether the window clips into the read from the left or right.

struct WindowClippingInfo_
{
    // Flag clipLeft is true if left window border clips into the read alignment.
    bool clipLeft;
    // Flag clipRight is true if right window border clips into the read alignment.
    bool clipRight;

    // Begin and end position in the read sequence for clipping.
    unsigned readBeginPos, readEndPos;
    // Begin and end position in the read alignment.
    unsigned readAliBeginPos, readAliEndPos;

    WindowClippingInfo_() : clipLeft(false), clipRight(false), readBeginPos(0), readEndPos(0),
                            readAliBeginPos(0), readAliEndPos(0)
    {}

    // Return true if window clipped into both borders.
    bool clipBoth() const
    {
        return clipLeft && clipRight;
    }
};


// ----------------------------------------------------------------------------
// Class AnsonMyersRealignmentRound_
// ----------------------------------------------------------------------------

// One iteration of the realignment.

template <typename TFragmentStore>
class AnsonMyersRealignmentRound_
{
public:
    typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
    typedef typename Value<TAlignedReadStore>::Type TAlignedReadElement;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type TAlignedReadIter;

    typedef typename TFragmentStore::TReadSeqStore      TReadSeqStore;
    typedef typename Value<TReadSeqStore>::Type         TReadSeq;
    typedef typename Value<TReadSeq>::Type              TStoreAlphabet;
    typedef typename BaseAlphabet<TStoreAlphabet>::Type TAlphabet;  // TODO(holtgrew): Rename?
    typedef ProfileChar<TAlphabet> TProfileChar;
    typedef String<TProfileChar>   TProfileString;
    typedef typename Iterator<TProfileString>::Type TProfileIter;

    // The FragmentStore for the sequences.
    TFragmentStore & store;
    // The contig's profile sequence.
    TProfileString & contigProfile;
    // The aligned reads for this contig.
    TAlignedReadStore & contigAlignedReads;
    // The timing helper struct.
    AnsonMyersTimes_ & times;
    // The options.
    RealignmentOptions_ const & options;

    // -----------------------------------------------------------------------
    // Public Interface
    // -----------------------------------------------------------------------

    AnsonMyersRealignmentRound_(TFragmentStore & store,
                                TProfileString & contigProfile,
                                TAlignedReadStore & contigAlignedReads,
                                AnsonMyersTimes_ & times,
                                RealignmentOptions_ const & options) :
            store(store), contigProfile(contigProfile), contigAlignedReads(contigAlignedReads), times(times),
            options(options)
    {}

    // Run one round of realignment.
    void run(unsigned windowBegin, unsigned & windowEnd);

    // Helper that returns whether the given profile character is for an all-gaps column.
    bool _allGaps(TProfileChar const & c)
    {
        return (c.count[0] + c.count[1] + c.count[2] + c.count[3] + c.count[4] == 0u);
    }

    bool _empty(TProfileChar const & c)
    {
        return (c.count[0] + c.count[1] + c.count[2] + c.count[3] + c.count[4] + c.count[5] == 0u);
    }

    // -----------------------------------------------------------------------
    // Implementation Functions
    // -----------------------------------------------------------------------

    inline bool _overlap(unsigned begin1, unsigned end1, unsigned begin2, unsigned end2)
    {
        return (begin1 < end2) && (begin2 < end1);
    }

    void _getReadAsProfile(TProfileString & profileStr,
                           TAlignedReadElement /*const*/ & el)
    {
        clear(profileStr);
        resize(profileStr, length(store.readSeqStore[el.readId]), TProfileChar());

        TProfileIter itP = begin(profileStr);

        typedef typename Iterator<TReadSeq, Standard>::Type TReadSeqIter;
        for (TReadSeqIter it = begin(store.readSeqStore[el.readId]),
                     itEnd = end(store.readSeqStore[el.readId]); it != itEnd; ++it, ++itP)
            itP->count[ordValue(*it)] = 1;
    }

    // Returns the information about the infix from the profile to extract.
    ExtractProfileInfo_ _getExtractPositions(TAlignedReadElement const & el)
    {
        ExtractProfileInfo_ result;
        if ((unsigned)el.beginPos < options.environment)
        {
            result.profBeginPos = 0;
            result.aliBeginPos = el.beginPos;
        }
        else
        {
            result.profBeginPos = el.beginPos - options.environment;
            result.aliBeginPos = options.environment;
        }

        if ((int)(el.endPos + options.environment) > (int)length(contigProfile))
            result.profEndPos = length(contigProfile);
        else
            result.profEndPos = el.endPos + options.environment;
        result.aliEndPos = result.aliBeginPos + (el.endPos - el.beginPos);

        // Sanity check on the result.
        SEQAN_ASSERT_LEQ(result.aliEndPos, result.profEndPos - result.profBeginPos);
        SEQAN_ASSERT_EQ((unsigned)(el.endPos - el.beginPos), result.aliEndPos - result.aliBeginPos);
        SEQAN_ASSERT_GEQ(result.profEndPos - result.profBeginPos, (unsigned)(el.endPos - el.beginPos));
        SEQAN_ASSERT_LEQ(result.profEndPos - result.profBeginPos,
                         (unsigned)(el.endPos - el.beginPos + 2 * options.environment));

        return result;
    }

    // Adjusts the ExtractProfileInfo_ object for the given window.
    void _updateExtractPositions(ExtractProfileInfo_ & info,
                                 unsigned windowBegin, unsigned windowEnd)
    {
        if (windowBegin == windowEnd)
            return;  // No update if window is disabled.

        // Get positions of read alignment from ExtractProfileInfo_.
        unsigned globalAliBeginPos = info.profBeginPos + info.aliBeginPos;
        unsigned globalAliEndPos = info.profBeginPos + info.aliEndPos;

        // Adjust begin positions in ExtractProfileInfo_.
        if (info.profBeginPos < windowBegin)
            info.profBeginPos = windowBegin;
        info.aliBeginPos = (globalAliBeginPos > info.profBeginPos) ? (globalAliBeginPos - info.profBeginPos) : 0;

        // Adjust end positions in ExtractProfileInfo_.
        if (info.profEndPos > windowEnd)
            info.profEndPos = windowEnd;
        if (globalAliEndPos > info.profEndPos)
            info.aliEndPos = info.profEndPos - info.profBeginPos;
        else
            info.aliEndPos = (info.profEndPos - info.profBeginPos) - (info.profEndPos - globalAliEndPos);
    }

    // Compute the WindowClippingInfo_.
    void _computeReadClippingInfo(WindowClippingInfo_ & windowInfo, TAlignedReadElement const & el,
                                  unsigned windowBegin, unsigned windowEnd)
    {
        if (windowBegin == windowEnd)  // no window
        {
            windowInfo.clipLeft = false;
            windowInfo.clipRight = false;
            windowInfo.readAliBeginPos = 0;
            windowInfo.readAliEndPos = el.endPos - el.beginPos;
            windowInfo.readBeginPos = 0;
            windowInfo.readEndPos = length(store.readSeqStore[el.readId]);
            return;
        }

        // Compute clipping flags.
        windowInfo.clipLeft = ((unsigned)el.beginPos < windowBegin);
        windowInfo.clipRight = ((unsigned)el.endPos > windowEnd);

        // Compute read alignment begin and end position.
        typedef Gaps<TReadSeq, AnchorGaps<String<typename TFragmentStore::TReadGapAnchor> > > TReadGaps;
        TReadGaps readGaps(store.readSeqStore[el.readId], el.gaps);
        if (options.debug)
            std::cerr << "READ GAPS\t" << readGaps << "\n";
        if (windowBegin < (unsigned)el.beginPos)
            windowInfo.readAliBeginPos = 0;
        else
            windowInfo.readAliBeginPos = windowBegin - el.beginPos;
        if (windowEnd > (unsigned)el.endPos)
            windowInfo.readAliEndPos = length(readGaps);
        else
            windowInfo.readAliEndPos = length(readGaps) - (el.endPos - windowEnd);

        // Compute source begin and end position.
        windowInfo.readBeginPos = toSourcePosition(readGaps, windowInfo.readAliBeginPos);
        windowInfo.readEndPos = toSourcePosition(readGaps, windowInfo.readAliEndPos);
    }

    // TODO(holtgrew): Build profile here to circumvent profileStr2 and the swap().
    // Subtract the counts for the alignment from the profile part.
    //
    // The aligned read element is given by copy since the original one can be altered through the removeGaps() call.
    void _subtractReadAlignment(unsigned windowBegin,
                                unsigned & windowEnd,
                                TProfileString & profileStr,
                                ExtractProfileInfo_ & info,
                                TAlignedReadElement el,
                                int elPos /*pos of el in aligned reads*/,
                                WindowClippingInfo_ const & windowInfo)
    {
        typedef typename Iterator<TProfileString, Standard>::Type TProfileIt;
        typedef Gaps<TReadSeq, AnchorGaps<String<typename TFragmentStore::TReadGapAnchor> > > TReadGaps;
        typedef typename Iterator<TReadGaps, Standard>::Type TReadGapsIter;

        TReadGaps readGaps(store.readSeqStore[el.readId], el.gaps);
        TReadGapsIter itR = iter(readGaps, windowInfo.readAliBeginPos, Standard());
        TReadGapsIter itREnd = iter(readGaps, windowInfo.readAliEndPos, Standard());
        TProfileIt itP = iter(profileStr, info.aliBeginPos, Standard());
        TProfileIt itPEnd = iter(profileStr, info.aliEndPos, Standard());
        (void)itPEnd;  // TODO(holtgrew): Still required?

        if (options.debug)
        {
            std::cerr << "READ GAPS\t";
            for (TReadGapsIter itR2 = itR; itR2 != itREnd; ++itR2)
                if (isGap(itR2))
                    std::cerr << "-";
                else
                    std::cerr << convert<char>(*itR2);
            std::cerr << "\n";

            std::cerr << "  ==> el.beginPos == " << el.beginPos << "\n";
        }

        // We will build a profileStr2 that has fewer gaps.
        TProfileString profileStr2;
        reserve(profileStr2, length(profileStr));
        // Append profile charactes without alignment by read.
        append(profileStr2, prefix(profileStr, info.aliBeginPos));

        // Subtract read alignment from profile, yielding profileStr2.
        int removeLeft = 0, removeRight = 0;  // number of bases removed left and right of read
        unsigned pos = info.aliBeginPos;   // position in profileStr
        unsigned posA = el.beginPos + windowInfo.readAliBeginPos;  // position in current multi-read aligment, including
                                                                   // removeGap()s.
        int gapsRemoved = 0;
        unsigned dbgCounter = 0;
        (void)dbgCounter;
        for (; itR != itREnd; ++itR, ++itP, ++pos)
        {
            SEQAN_ASSERT_MSG(itP != end(profileStr, Standard()),
                             "Profile part includes read alignment (debug counter = %u).", dbgCounter);

            unsigned idx = isGap(itR) ? valueSize<TAlphabet>() : ordValue(*itR);
#if SEQAN_ENABLE_DEBUG
            int posP = itP - begin(profileStr, Standard());
            int posR = itR - begin(readGaps, Standard());
#endif  // #if SEQAN_ENABLE_DEBUG
            // std::cerr << "idx == " << idx << "\tpos == " << pos << "\tidx == " << idx << "\tposP == " << posP
            //           << "\tposR == " << posR << "\t" << "posA == " << posA << "\t"
            //           << "itP->count[0] == " << itP->count[0] << "\t" << "itP->count[1] == " << itP->count[1] << "\t"
            //           << "itP->count[2] == " << itP->count[2] << "\t" << "itP->count[3] == " << itP->count[3] << "\t"
            //           << "itP->count[4] == " << itP->count[4] << "\t" << "itP->count[5] == " << itP->count[5] << "\t"
            //           << "*itR == " << (isGap(itR) ? '-' : convert<char>(TAlphabet(*itR))) << "\n";
            SEQAN_ASSERT_GT_MSG(itP->count[idx], 0u,
                                "Must have count for alignment char (pos = %u, idx = %u, posP = %d, posR = %d).",
                                pos, idx, posP, posR);
            itP->count[idx] -= 1;

            // Append to profile2 unless there is an all-gaps column and there is a gap in the read gaps.
            if ((!_allGaps(*itP) || isGap(itR)) && !_empty(*itP))
            {
                appendValue(profileStr2, *itP);
                ++posA;
            }
            else
            {
                if (options.debug)
                {
                    swap(store.alignedReadStore, contigAlignedReads);
                    std::cerr << "ALIGNMENT is (current == " << store.readSeqStore[el.readId] << ")\n";
                    seqan::AlignedReadLayout layout;
                    layoutAlignment(layout, store);
                    std::stringstream ss;
                    printAlignment(std::cerr, layout, store, front(store.alignedReadStore).contigId,
                                   0, 1000, 0, 100);
                    swap(store.alignedReadStore, contigAlignedReads);

                    std::cerr << " => removeGap(contigAlignedReads, " << posA << ", " << elPos << ")\n"
                              << "       pos == " << pos << "\n";
                    // std::cerr << "    itP:";
                    // for (unsigned i = 0; i < 6u; ++i)
                    //     std::cerr << "\titP->count[" << i << "] == " << itP->count[i];
                    // std::cerr << "\n";
                }
                _removeGap(contigAlignedReads, posA, iter(contigAlignedReads, elPos, Standard()));
                gapsRemoved += 1;  // WAS: (removeGap(contigAlignedReads, posA, iter(contigAlignedReads, elPos,
                                   //                 Standard())) > 0);
                if (pos < info.aliBeginPos)
                    removeLeft += 1;  // TODO(holtgrew): cannot be reached!
                if (pos >= info.aliEndPos)
                    removeRight += 1;  // TODO(holtgrew): can also not be reached!
            }
        }

        if (options.debug)
            std::cerr << "gapsRemoved == " << gapsRemoved << "\n"
                      << "adjusting window from " << windowBegin << " " << windowEnd << "\n";
        if (windowBegin != windowEnd)
            windowEnd -= gapsRemoved;
        if (options.debug)
            std::cerr << "adjusting window to " << windowBegin << " " << windowEnd << "\n";

        // Add remaining profile characters.
        append(profileStr2, suffix(profileStr, pos));

        SEQAN_ASSERT_LEQ(removeLeft, (int)info.aliBeginPos);
        SEQAN_ASSERT_LEQ(removeRight, (int)(length(profileStr) - info.aliEndPos));
        info.aliBeginPos -= removeLeft;
        info.aliEndPos -= removeRight;

        // Write out other profile string.
        swap(profileStr, profileStr2);
    }

    // Use the alignment result in matches to update contigProfile and insert/remove gaps from the fragment store.
    //
    // matches     -- alignment result
    // profilePart -- the part of the profile that was extracted
    // info        -- extraction information about profilePart
    // el          -- aligned read store element
    void _updateAlignments(
            unsigned windowBegin,
            unsigned & windowEnd,  // only updated
            TProfileString & profilePart,
            // Align<TProfileString> /*const*/ & align,
            Gaps<TProfileString> /*const*/ & profileGaps,
            Gaps<TReadSeq> /*const*/ & readGaps,
            ExtractProfileInfo_ & info,
            WindowClippingInfo_ const & windowInfo,
            TAlignedReadIter const & it);
};

// ----------------------------------------------------------------------------
// Class AnsonMyersRealigner_
// ----------------------------------------------------------------------------

// Implementation of the Anson-Myers algorithm for realigning the fragments in a FragmentStore.
//
// This class is used for holding references to the input data structures and the helper data structures.

template <typename TFragmentStore>
class AnsonMyersRealigner_
{
public:
    typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type TAlignedReadIter;

    typedef typename TFragmentStore::TReadSeq           TReadSeq;
    typedef typename Value<TReadSeq>::Type              TStoreAlphabet;
    typedef typename BaseAlphabet<TStoreAlphabet>::Type TAlphabet;  // TODO(holtgrew): Rename?
    typedef ProfileChar<TAlphabet> TProfileChar;
    typedef String<TProfileChar>   TProfileString;
    typedef typename Iterator<TProfileString, Standard>::Type TProfileIter;

    // The FragmentStore object with the reads.
    TFragmentStore & store;
    // The configuration to use for the realignment.
    RealignmentOptions_ options;

    // A copy of the aligned reads for the current contig.
    TAlignedReadStore contigAlignedReads;
    // Information on the alignments of the current contig.
    ContigAlignmentsInfo_<TFragmentStore> contigAlignmentInfos;

    // The contig's profile sequence.
    TProfileString contigProfile;

    // Timing information.
    AnsonMyersTimes_ times;

    // -----------------------------------------------------------------------
    // Public Interface
    // -----------------------------------------------------------------------

    // Construct
    AnsonMyersRealigner_(TFragmentStore & store, RealignmentOptions_ const & options) :
            store(store), options(options)
    {}

    // Run the realignment on the given contig.  If windowBegin != windowEnd then only the given window is realigned.
    // The value of windowEnd is adjusted in case columns are added to or removed from the window.
    void run(unsigned contigID, unsigned windowBegin, unsigned &windowEnd);

    // Run on the whole contig without a window.
    void run(unsigned contigID)
    {
        unsigned windowBegin = 0;
        unsigned windowEnd = 0;
        run(contigID, windowBegin, windowEnd);
    }

    // -----------------------------------------------------------------------
    // Implementation Functions
    // -----------------------------------------------------------------------

    // Check read alignments for consistency.
    void _checkReadAlignments(TAlignedReadStore & alignedReadStore)
    {
        if (options.debug)
            std::cerr << "begin/end positions of alignments (central check)\n";
        TAlignedReadIter itEnd = end(alignedReadStore, Standard());
        for (TAlignedReadIter it2 = begin(alignedReadStore, Standard()); it2 != itEnd; ++it2)
        {
            typedef Gaps<TReadSeq, AnchorGaps<String<typename TFragmentStore::TReadGapAnchor> > > TReadGaps;
            TReadGaps readGaps(store.readSeqStore[it2->readId], it2->gaps);
            // if (options.debug)
            //     std::cerr << "it2->beginPos == " << it2->beginPos << ", it2->endPos == " << it2->endPos << "\n";
            if ((unsigned)abs((int)(it2->endPos - it2->beginPos)) != length(readGaps))
            {
                std::cerr << "READ GAPS\t>>>" << readGaps << "<<<\n";
                SEQAN_FAIL("Inconsistent begin/endPos");
            }
        }
    }

    void _checkReadAlignments()
    {
        _checkReadAlignments(store.alignedReadStore);
        _checkReadAlignments(contigAlignedReads);
    }

    // Sort aligned reads, fill contigAlignedReads and fill contigAlignmentInfos, add reference sequence as read if
    // configured to do so.
    void _beginContig(unsigned contigID);

    // Write back the updated information from contigAlignedReads to store, add reference sequence alignment if
    // configured to do so.
    void _endContig(unsigned contigID);
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function AnsonMyersRealigner_::run()
// ----------------------------------------------------------------------------

template <typename TProfileGaps, typename TReadGaps, typename TConsensusScore, typename TAlignConfig,
          typename TAlgoTag>
void _readToProfileAlignment(TProfileGaps & profileGaps, TReadGaps & readGaps, TConsensusScore const & consScore,
                             TAlignConfig const & alignConfig, int & lowerBand, int & upperBand, TAlgoTag const & tag)
{
    _fixBandSize(lowerBand, upperBand, source(profileGaps), source(readGaps), alignConfig, tag);
    globalAlignment(profileGaps, readGaps, consScore, alignConfig, /*lowerBand, upperBand, */tag);
}

template <typename TProfileGaps, typename TReadGaps, typename TConsensusScore, typename TAlignConfig>
void _readToProfileAlignment(TProfileGaps & profileGaps, TReadGaps & readGaps, TConsensusScore const & consScore,
                             TAlignConfig const & alignConfig, int & lowerBand, int & upperBand, bool linear)
{
    // Handle the case where the read was the only one aligning at this position.  It would be impossible to compute a
    // pairwise alignment below.
    if (empty(source(profileGaps)))
    {
        clearGaps(readGaps);
        clearGaps(profileGaps);
        insertGaps(profileGaps, 0, length(readGaps));
    }
    else
    {
        if (linear)
            _readToProfileAlignment(profileGaps, readGaps, consScore, alignConfig, lowerBand, upperBand,
                                    NeedlemanWunsch());
        else
            _readToProfileAlignment(profileGaps, readGaps, consScore, alignConfig, lowerBand, upperBand,
                                    Gotoh());
    }
}

template <typename TFragmentStore>
void AnsonMyersRealigner_<TFragmentStore>::run(unsigned contigID, unsigned windowBegin, unsigned & windowEnd)
{
    _checkReadAlignments();

    // Fill contigAlignedReads, contigAlignmentInfos, contigProfile.
    _beginContig(contigID);

    _checkReadAlignments();

    // Perform realignment round until the score does not decrease any more.
    int score = maxValue<int>();
    int oldScore = maxValue<int>();
    for (unsigned roundNo = 0; score == maxValue<int>() || score < oldScore; ++roundNo)
    {
        if (options.debug)
            std::cerr << "Realignment round " << roundNo << "\n";
        oldScore = score;
        AnsonMyersRealignmentRound_<TFragmentStore> round(store, contigProfile, contigAlignedReads, times, options);
        _checkReadAlignments();
        round.run(windowBegin, windowEnd);
        _checkReadAlignments();
        // TODO(holtgrew): This is not a similarity score but more of a distance score.
        score = scoreConsensus(contigProfile);
        if (options.debug)
            std::cerr << "Old score = " << oldScore << ", score == " << score << "\n";
    }

    _checkReadAlignments();

    // Write back the realigned read alignments into store and set contig to the consensus in case that the original
    // reference was included as a pseudo-read.
    _endContig(contigID);

    _checkReadAlignments();
}

// ----------------------------------------------------------------------------
// Function AnsonMyersRealigner_::_beginContig()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Completely restore RC state of sequence and begin/end position again.
// TODO(holtgrew): Split into smaller functions.

template <typename TFragmentStore>
void AnsonMyersRealigner_<TFragmentStore>::_beginContig(unsigned contigID)
{
    double startTime = sysTime();
    contigAlignmentInfos.reset();  // reset/clear

    // ------------------------------------------------------------------------
    // Read Sorting
    // ------------------------------------------------------------------------

    // Sort aligned reads by contig id and get iterators to begin/end of the alignments on the contig.
    sortAlignedReads(store.alignedReadStore, SortContigId());
    contigAlignmentInfos.alignedItBegin = lowerBoundAlignedReads(store.alignedReadStore, contigID, SortContigId());
    contigAlignmentInfos.alignedItEnd = upperBoundAlignedReads(store.alignedReadStore, contigID, SortContigId());

    // Sort the reads on the current contig according to their begin posiiton.
    sortAlignedReads(infix(store.alignedReadStore,
                           contigAlignmentInfos.alignedItBegin - begin(store.alignedReadStore, Standard()),
                           contigAlignmentInfos.alignedItEnd - begin(store.alignedReadStore, Standard())),
                     SortBeginPos());

    // ------------------------------------------------------------------------
    // Extract Alignment; Reverse-Complement Fragment Sequence
    // ------------------------------------------------------------------------

    // Copy out the read alignments belonging to the contig and reverse-complement the read sequence if necessary.
    //
    // We will later work on the copy of the read alignments and at the end we will write it back to store.  The
    // reverse-complements are also used in the alignments.  In _endContig(), we will restore the original sequence
    // again.
    //
    // TODO(holtgrew): We could also work on an infix of the read alignments. Advantages/Disadvantages?
    clear(contigAlignedReads);
    typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
    typedef typename Value<TAlignedReadStore>::Type TAlignedElement;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type TAlignedReadIter;
    for (TAlignedReadIter it = contigAlignmentInfos.alignedItBegin; it != contigAlignmentInfos.alignedItEnd; ++it)
    {
        TAlignedElement el = *it;
        if (it->beginPos > it->endPos)  // reverse-complemented
        {
            reverseComplement(store.readSeqStore[it->readId]);
            std::swap(el.beginPos, el.endPos);
        }
        contigAlignmentInfos.minPos = std::min(contigAlignmentInfos.minPos, el.beginPos);
        contigAlignmentInfos.maxPos = std::max(contigAlignmentInfos.maxPos, el.endPos);
        appendValue(contigAlignedReads, el);
    }

    // ------------------------------------------------------------------------
    // Include Reference Sequence
    // ------------------------------------------------------------------------

    if (options.includeReference)
    {
        // Compute read and match ID.
        unsigned readID = length(store.readSeqStore);
        unsigned matchID = length(store.alignedReadStore);

        // Append sequence and name of pseudo-read.
        appendRead(store, store.contigStore[contigID].seq);
        appendValue(store.readNameStore, store.contigNameStore[contigID]);
        store.contigNameStore[contigID] += "Consensus_";  // TODO(holtgrew): Prepend?

        // Create alignment for pseudo-read.
        typedef typename TFragmentStore::TReadPos       TReadPos;
        typedef typename TFragmentStore::TContigStore   TContigStore;
        typedef typename Value<TContigStore>::Type      TContig;
        // typedef typename TFragmentStore::TContigPos     TContigPos;
        typedef typename TFragmentStore::TContigSeq     TContigSeq;
        typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> > TContigGaps;
        // Construct aligned read element.
        TAlignedElement el;
        el.id = matchID;
        el.readId = readID;
        el.contigId = contigID;
        contigAlignmentInfos.minPos = el.beginPos = 0;
        el.gaps = store.contigStore[contigID].gaps;
        // Compute begin/end positions and append element to the contig's read alignments.
        TContigGaps contigGaps(store.contigStore[contigID].seq, store.contigStore[contigID].gaps);
        TReadPos contigLength = length(store.contigStore[contigID].seq);
        TReadPos contigEndGapPos = positionSeqToGap(contigGaps, contigLength - 1) + 1;
        contigAlignmentInfos.maxPos = el.endPos = std::max(contigAlignmentInfos.maxPos,
                                                           std::max(contigEndGapPos, contigLength));
        appendValue(contigAlignedReads, el);
    }

    // ------------------------------------------------------------------------
    // Create Profile Sequence
    // ------------------------------------------------------------------------

    clear(contigProfile);
    resize(contigProfile, contigAlignmentInfos.maxPos - contigAlignmentInfos.minPos, TProfileChar());

    // Iterate over the pairwise alignments of the reads against the contig that are stored in the fragment store.
    for (TAlignedReadIter it = begin(contigAlignedReads, Standard());
         it != end(contigAlignedReads, Standard()); ++it)
    {
        // Shift the aligned reads from the contig to the left by contigAlignmentInfos.minPos;
        it->beginPos -= contigAlignmentInfos.minPos;
        it->endPos -= contigAlignmentInfos.minPos;

        typedef Gaps<TReadSeq, AnchorGaps<String<typename TFragmentStore::TReadGapAnchor> > > TReadGaps;
        typedef typename Iterator<TReadGaps, Standard>::Type TReadGapsIter;

        // We create Gaps for the read we are iterating over.
        // TODO(holtgrew): Is clipping stored in the read gaps? It appears so in Tobias' original code.
        TReadGaps readGaps(store.readSeqStore[it->readId], it->gaps);

        // Get iterators to beginning of read gaps and the corresponding positions in the contig gaps and the profile
        // string.
        TProfileIter itP = iter(contigProfile, it->beginPos, Standard());
        TReadGapsIter itR = begin(readGaps, Standard());
        TReadGapsIter itREnd = end(readGaps, Standard());

        // TODO(holtgrew): This could become easier if there was a countTrailingGaps() function.
        // Skip over leading gaps, advance iterator into profile.
        itP += countGaps(itR);
        itR += countGaps(itR);
        // Skip over trailing gaps.
        if (itR != itREnd)
        {
            --itREnd;
            while (itREnd != itR && isGap(itREnd))
                --itREnd;
            if (!isGap(itREnd))
                ++itREnd;
        }

        // Copy over characters/gaps into profile.
        for (; itR != itREnd; ++itR, ++itP)
            ++itP->count[isGap(itR) ? valueSize<TAlphabet>() : ordValue(*itR)];
    }

    if (options.debug)
    {
        std::cerr << "INITIAL PROFILE\n";
        _printProfile(std::cerr, contigProfile);
    }

    // Update times used for processing beginning of contig.
    times.mainBeginContig += sysTime() - startTime;
}

// ----------------------------------------------------------------------------
// Function AnsonMyersRealigner_::_endContig()
// ----------------------------------------------------------------------------

template <typename TFragmentStore>
void AnsonMyersRealigner_<TFragmentStore>::_endContig(unsigned contigID)
{
    // -----------------------------------------------------------------------
    // Update Read Alignments
    // -----------------------------------------------------------------------

    TAlignedReadIter it = contigAlignmentInfos.alignedItBegin;
    TAlignedReadIter itEnd = contigAlignmentInfos.alignedItEnd;
    TAlignedReadIter itNew = begin(contigAlignedReads, Standard());

    for (unsigned idx = 0; it != itEnd; ++idx, ++it, ++itNew)
    {
        it->beginPos = itNew->beginPos;
        it->endPos = itNew->endPos;
        it->gaps = itNew->gaps;  // TODO(holtgrew): Old code removed empty anchors -- required?
        if (it->endPos < it->beginPos)
            reverseComplement(store.readSeqStore[it->readId]);
    }

    // -----------------------------------------------------------------------
    // Update Contig
    // -----------------------------------------------------------------------

    // The contig will be replaced by the consensus sequence.

    typedef typename TFragmentStore::TContigStore   TContigStore;
    typedef typename Value<TContigStore>::Type      TContig;
    typedef typename TContig::TContigSeq            TContigSeq;
    typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> > TContigGaps;
    typedef typename Position<TContigSeq>::Type     TContigPos;

    // typedef typename Iterator<TContigGaps, Standard>::Type TContigGapsIter;

    // Get shortcut to contig sequence and clear for rebuilding.
    TContigSeq & contigSeq = store.contigStore[contigID].seq;
    clear(contigSeq);

    // Iterate over the profile and build up the contig sequence and store gap positions.
    unsigned gapLen = 0;  // length of current gap
    TContigPos contigPos = 0;  // with gaps
    String<Pair<unsigned, unsigned> > gaps;  // Gap positions and lengths.
    for (TProfileIter it = begin(contigProfile, Standard()), itEnd = end(contigProfile, Standard());
         it != itEnd; ++it, ++contigPos)
    {
        if (static_cast<char>(*it) == gapValue<char>())
        {
            ++gapLen;
            continue;
        }

        // Record gaps position.
        if (gapLen)
            appendValue(gaps, Pair<unsigned, unsigned>(contigPos - gapLen, gapLen));
        gapLen = 0;

        // Append to contig sequence.
        //
        // TODO(weese): Here we convert from ProfileChar<Dna5>->Dna5->Dna5Q
        // instead diverting through Dna5 we could think of directly converting
        // a profile to a quality value, e.g. like the base caller Phred does.
        // Therefore the conversion ProfileChar<Dna5> <-> Dna5Q needs to be
        // defined.
        appendValue(contigSeq, static_cast<TAlphabet>(*it));
    }

    // Rebuild contig gaps.
    TContigGaps contigGaps(store.contigStore[contigID].seq, store.contigStore[contigID].gaps);
    clearGaps(contigGaps);
    for (unsigned i = 0; i < length(gaps); ++i)
    {
        if (options.debug)
            std::cerr << "insertGaps(contigGaps, " << gaps[i].i1 << ", " << gaps[i].i2 << ")\n";
        insertGaps(contigGaps, gaps[i].i1, gaps[i].i2);
    }

    if (options.debug)
    {
        std::cerr << "Resulting contig gaps\n"
                  << contigGaps << "\n";
    }

    // -----------------------------------------------------------------------
    // Include Reference
    // -----------------------------------------------------------------------

    if (options.includeReference)
        appendValue(store.alignedReadStore, back(contigAlignedReads));
}


// ----------------------------------------------------------------------------
// Function _fixBandSize()
// ----------------------------------------------------------------------------

/*
 * @fn _fixBandSize
 * @headerfile <seqan/align.h>
 * @brief Fix a band for alignment given sequences (for their lengths) and an @link AlignConfig @endlink.
 *
 * @signature void _fixBandSize(lowerDiag, upperDiag, seqH, seqV, alignConfig);
 *
 * @param[in,out] lowerDiag   The lower band position (<tt>int</tt>) to adjust.
 * @param[in,out] upperDiag   The upper band position (<tt>int</tt>) to adjust.
 * @param[in]     seqH        The @link ContainerConcept container @endlink to use in horizontal direction of alignment
 *                            matrix.
 * @param[in]     seqV        The @link ContainerConcept container @endlink to use in vertical direction of alignment
 *                            matrix.
 * @param[in]     alignConfig The @link AlignConfig @endlink to use for the alignment configuration.
 */

template <typename TSequenceH, typename TSequenceV, typename TAlignConfig, typename TAlgoTag>
void _fixBandSize(int & lDiag,
                  int & uDiag,
                  TSequenceH const & seqH,
                  TSequenceV const & seqV,
                  TAlignConfig const & /*alignConfig*/,
                  TAlgoTag const & /*algoTag*/)
{
    typedef typename SubstituteAlignConfig_<TAlignConfig>::Type TFreeEndGaps;
    typedef typename If<typename IsSameType<TAlgoTag, Gotoh>::Type, AffineGaps, LinearGaps>::Type TGapsType;
    typedef typename SetupAlignmentProfile_<DPGlobal, TFreeEndGaps, TGapsType,
                                            TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > >::Type TDPProfile;

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

// ----------------------------------------------------------------------------
// Function AnsonMyersRealignmentRound_::run()
// ----------------------------------------------------------------------------

template <typename TFragmentStore>
void AnsonMyersRealignmentRound_<TFragmentStore>::run(unsigned windowBegin, unsigned & windowEnd)
{
    if (empty(contigAlignedReads))
        return;
    // If there is only one read on this contig and no window is set or the window spans the whole alignment then remove
    // all gaps from this one read and exit.
    if (length(contigAlignedReads) == 1u)  // only one alignment
    {
        if (((windowBegin == 0u) && (windowBegin == windowEnd)) ||  // no window
            (windowBegin <= (unsigned)contigAlignedReads[0].beginPos &&
             windowEnd >= (unsigned)contigAlignedReads[0].endPos))  // spans only
                                                                                                           // alignment
        {
            contigAlignedReads[0].beginPos = 0;
            contigAlignedReads[0].endPos = length(store.readSeqStore[contigAlignedReads[0].readId]);
            clear(contigAlignedReads[0].gaps);
            _getReadAsProfile(contigProfile, contigAlignedReads[0]);
            return;
        }
    }

    double roundStartTime = sysTime();  // round start time
    if (options.debug)
        std::cerr << "STARTING ROUND WITH windowBegin == " << windowBegin << ", windowEnd == " << windowEnd << "\n";
    SEQAN_ASSERT_LEQ(windowBegin, windowEnd);
    if (windowBegin == windowEnd && windowBegin != 0u)
    {
        if (options.debug)
            std::cerr << "Window is empty.\n";
        return;
    }

    // The currently extracted part from the sequence and the read as profile (ord values stored in read[i].count[0].
    TProfileString profilePart;

    // Remove each read from the profile, realign it to the profile, and update the profile.
    TAlignedReadIter itEnd = end(contigAlignedReads, Standard());
    for (TAlignedReadIter it = begin(contigAlignedReads, Standard()); it != itEnd; ++it)
    {
        if (options.includeReference && it->readId + 1 == length(store.readSeqStore))
            continue;  // do not subtract reference pseudo-read
        // Ignore reads that do not overlap with the window.
        if (windowBegin != windowEnd && !_overlap(windowBegin, windowEnd, it->beginPos, it->endPos))
            continue;

        // Extract the profile around the alignment.
        double beginExtractProfile = sysTime();
        // Get positions in the alignment to extrad and positions in the alignment in this part.
        ExtractProfileInfo_ info = _getExtractPositions(*it);
        SEQAN_ASSERT_EQ(info.profBeginPos + info.aliBeginPos, (unsigned)it->beginPos);
        SEQAN_ASSERT_EQ(info.profBeginPos + info.aliEndPos, (unsigned)it->endPos);
        if (options.debug)
            std::cerr << "WINDOW                       [" << windowBegin << ", " << windowEnd << ")\n"
                      << "UNADJUSTED\n"
                      << "Extracting profile begin/end " << info.profBeginPos << "/" << info.profEndPos << "\n"
                      << "Aligned read begin/end       " << it->beginPos << "/" << it->endPos << "\n"
                      << "Alignment begin/end          " << info.aliBeginPos << "/" << info.aliEndPos << "\n"
                      << "Read sequence                " << store.readSeqStore[it->readId] << "\n"
                      << "it->beginPos == " << it->beginPos << ", it->endPos == " << it->endPos << "\n";
        // Update these extraction positions for the given window.
        _updateExtractPositions(info, windowBegin, windowEnd);
        // Compute read window clipping info.
        WindowClippingInfo_ windowInfo;
        _computeReadClippingInfo(windowInfo, *it, windowBegin, windowEnd);
        if (options.debug)
            std::cerr << "ADJUSTED\n"
                      << "Extracting profile begin/end " << info.profBeginPos << "/" << info.profEndPos << "\n"
                      << "Aligned read begin/end       " << it->beginPos << "/" << it->endPos << "\n"
                      << "Alignment begin/end          " << info.aliBeginPos << "/" << info.aliEndPos << "\n"
                      << "Read sequence                " << store.readSeqStore[it->readId] << "\n"
                      << "window clips left?  " << windowInfo.clipLeft << "\n"
                      << "window clips right? " << windowInfo.clipRight << "\n"
                      << "read begin pos      " << windowInfo.readBeginPos << "\n"
                      << "read end pos        " << windowInfo.readEndPos << "\n"
                      << "read ali begin pos  " << windowInfo.readAliBeginPos << "\n"
                      << "read ali end pos    " << windowInfo.readAliEndPos << "\n";
        if (windowInfo.readBeginPos == windowInfo.readEndPos)
            continue;  // Skip if we only cut out gap characters.
        SEQAN_ASSERT_EQ(windowInfo.readAliEndPos - windowInfo.readAliBeginPos,
                        info.aliEndPos - info.aliBeginPos);
        // Extract part of the profile with buffer bases left and right.
        profilePart = infix(contigProfile, info.profBeginPos, info.profEndPos);
        if (options.debug)
        {
            typedef typename TFragmentStore::TContigStore   TContigStore;
            typedef typename Value<TContigStore>::Type      TContig;
            typedef typename TContig::TContigSeq            TContigSeq;
            typedef typename Position<TContigSeq>::Type TContigPos;
            TContigSeq contigSeq;

            unsigned gapLen = 0;  // length of current gap
            TContigPos contigPos = 0;  // with gaps
            String<Pair<unsigned, unsigned> > gaps;  // Gap positions and lengths.
            for (TProfileIter it2 = begin(contigProfile, Standard()), itEnd = end(contigProfile, Standard());
                 it2 != itEnd; ++it2, ++contigPos)
            {
                if (static_cast<char>(*it2) == gapValue<char>())
                {
                    ++gapLen;
                    continue;
                }

                // Record gaps position.
                if (gapLen)
                    appendValue(gaps, Pair<unsigned, unsigned>(contigPos - gapLen, gapLen));
                gapLen = 0;

                // Append to contig sequence.
                //
                // TODO(weese): Here we convert from ProfileChar<Dna5>->Dna5->Dna5Q
                // instead diverting through Dna5 we could think of directly converting
                // a profile to a quality value, e.g. like the base caller Phred does.
                // Therefore the conversion ProfileChar<Dna5> <-> Dna5Q needs to be
                // defined.
                appendValue(contigSeq, static_cast<TAlphabet>(*it2));
            }
            std::cerr << "Consensus                    " << contigSeq << "\n";

            swap(store.alignedReadStore, contigAlignedReads);
            std::cerr << "ALIGNMENT is\n";
            seqan::AlignedReadLayout layout;
            layoutAlignment(layout, store);
            std::stringstream ss;
            printAlignment(std::cerr, layout, store, front(store.alignedReadStore).contigId,
                           0, 1000, 0, 100);
            swap(store.alignedReadStore, contigAlignedReads);

            std::cerr << "PROFILE\n";
            _printProfile(std::cerr, contigProfile);
            std::cerr << "PROFILE PART BEFORE SUBTRACTION\n";
            _printProfile(std::cerr, profilePart);
        }

        // Subtract the read alignment from the profile.
        _subtractReadAlignment(windowBegin, windowEnd, profilePart, info, *it, it - begin(contigAlignedReads,
                                                                                          Standard()),
                               windowInfo);
        SEQAN_ASSERT_LEQ(windowBegin, windowEnd);
        times.realignExtractProfile += sysTime() - beginExtractProfile;

        if (options.debug)
        {
            std::cerr << "PROFILE PART AFTER SUBTRACTION\n";
            _printProfile(std::cerr, profilePart);

            swap(store.alignedReadStore, contigAlignedReads);
            std::cerr << "AFTER SUBTRACTION of " << store.readSeqStore[it->readId] << " id= " << it->readId << "\n";
            seqan::AlignedReadLayout layout;
            layoutAlignment(layout, store);
            std::stringstream ss;
            printAlignment(std::cerr, layout, store, front(store.alignedReadStore).contigId,
                           0, 1000, 0, 100);
            swap(store.alignedReadStore, contigAlignedReads);
        }

        // Realign the read against the profile.
        double beginAlign = sysTime();
        Gaps<TProfileString> profileGaps(profilePart);
        Gaps<TReadSeq> readGaps;
        assignSource(readGaps, infix(store.readSeqStore[it->readId], windowInfo.readBeginPos, windowInfo.readEndPos));
        if (options.debug)
        {
            std::cerr << "READ PART; SOURCE OF READ GAPS (seq = " << store.readSeqStore[it->readId] << ")\n";
            std::cerr << source(readGaps) << "\n";
        }

        // Setup the score (same as in Anson-Myers paper), requires knowledge about the profile that is to be aligned
        // against.
        Score<int, WeightedConsensusScore<Score<int, ProfileSeqFracScore>,
                                          Score<int, ProfileSeqScore> > > consScore;
        assignProfile(consScore, profilePart);

        // Compute the flag isSingleLeft that indicates whether the current read is the first read and the coverage at
        // the first position is 1 (i.e. only this reads aligns there).  In this case, we have to use an overlap
        // alignment here instead of a semiglobal alignment and also adjust the band for the alignment.
        bool isSingleLeft = (it->beginPos == 0u) &&
                (it == begin(contigAlignedReads, Standard())) &&
                (it - begin(contigAlignedReads, Standard()) + 1 < (int)length(contigAlignedReads)) &&
                (((it + 1)->beginPos - it->beginPos) > 0);
        if (windowBegin != windowEnd && windowBegin > 0)
            isSingleLeft = false;  // Cannot be case of single left if window does not start at begin.
        int bandDelta = isSingleLeft ? ((it + 1)->beginPos - it->beginPos) : 0;
        // Compute the upper and lower band for the alignment of the read to the profile.  The band will be fixed
        // further using _fixBandSize() below.
        int upperBand = std::max((int)info.aliBeginPos + (int)options.bandwidth + bandDelta,
                                 (int)length(profilePart) - (int)length(readGaps) + (int)options.bandwidth + bandDelta);
        int lowerBand = std::min(-(int)options.bandwidth - bandDelta,
                                 (int)length(profilePart) - (int)length(readGaps) - (int)options.bandwidth - bandDelta);

        if (options.debug)
            std::cerr << "BEFORE ALIGNMENT\n"
                      << "length of align row 0\t" << length(profileGaps) << "\n"
                      << "length of align row 1\t" << length(readGaps) << "\n"
                      << "\n"
                      << "PROFILE\t" << profileGaps << "\n"
                      << "READ   \t" << readGaps << "\n";

        // TODO(holtgrew): Adjust begin/end gap scoring.
        bool linear = (options.method == RealignmentOptions_::ANSON_MYERS_NW);
        if (isSingleLeft)
        {
            if (options.debug)
                std::cerr << "CASE OF LEFT (lowerBand == " << lowerBand << ", upperBand == " << upperBand
                          << ", bandDelta == " << bandDelta << ")\n";
            AlignConfig<false, true, false, true> alignConfig;
            _readToProfileAlignment(profileGaps, readGaps, consScore, alignConfig, lowerBand, upperBand, linear);
        }
        else
        {
            if (windowInfo.clipLeft)
            {
                if (options.debug)
                    std::cerr << "Clipped left. (lowerBand == " << lowerBand << ", upperBand == " << upperBand << ")\n";
                AlignConfig<false, false, false, true> alignConfig;
                _readToProfileAlignment(profileGaps, readGaps, consScore, alignConfig, lowerBand, upperBand, linear);
            }
            else if (windowInfo.clipRight)
            {
                if (options.debug)
                    std::cerr << "Clipped right. (lowerBand == " << lowerBand << ", upperBand == "
                              << upperBand << ")\n";
                AlignConfig<true, false, false, false> alignConfig;
                _readToProfileAlignment(profileGaps, readGaps, consScore, alignConfig, lowerBand, upperBand, linear);
            }
            else if (windowInfo.clipBoth())
            {
                if (options.debug)
                    std::cerr << "Clipped both. (lowerBand == " << lowerBand << ", upperBand == "
                              << upperBand << ")\n";
                AlignConfig<false, false, false, false> alignConfig;
                _readToProfileAlignment(profileGaps, readGaps, consScore, alignConfig, lowerBand, upperBand, linear);
            }
            else // clipNone
            {
                if (options.debug)
                    std::cerr << "NOT case of left (lowerBand == " << lowerBand << ", upperBand == "
                              << upperBand << ")\n";
                AlignConfig<true, false, false, true> alignConfig;
                _readToProfileAlignment(profileGaps, readGaps, consScore, alignConfig, lowerBand, upperBand, linear);
            }
        }
#if 0
        if (isSingleLeft)
        {
            AlignConfig<false, true, false, true> alignConfig;
            if (options.debug)
                std::cerr << "CASE OF LEFT (lowerBand == " << lowerBand << ", upperBand == " << upperBand
                          << ", bandDelta == " << bandDelta << ")\n";
            if (options.method == RealignmentOptions_::ANSON_MYERS_NW)
            {
                _fixBandSize(lowerBand, upperBand, profilePart, source(readGaps), alignConfig, NeedlemanWunsch());
                globalAlignment(profileGaps, readGaps, consScore, alignConfig, lowerBand, upperBand,
                                NeedlemanWunsch());
            }
            else
            {
                _fixBandSize(lowerBand, upperBand, profilePart, source(readGaps), alignConfig, Gotoh());
                globalAlignment(profileGaps, readGaps, consScore, alignConfig, lowerBand, upperBand, Gotoh());
            }
        }
        else
        {
            AlignConfig<true, false, false, true> alignConfig;
            if (options.debug)
                std::cerr << "NOT case of left (lowerBand == " << lowerBand << ", upperBand == "
                          << upperBand << ")\n";
            if (options.method == RealignmentOptions_::ANSON_MYERS_NW)
            {
                _fixBandSize(lowerBand, upperBand, profilePart, source(readGaps), alignConfig, NeedlemanWunsch());
                globalAlignment(profileGaps, readGaps, consScore, alignConfig, lowerBand, upperBand,
                                NeedlemanWunsch());
            }
            else
            {
                _fixBandSize(lowerBand, upperBand, profilePart, source(readGaps), alignConfig, Gotoh());
                globalAlignment(profileGaps, readGaps, consScore, alignConfig, lowerBand, upperBand, Gotoh());
            }
        }
#endif  // #if 0
        times.realignAlign = sysTime() - beginAlign;

        if (options.debug)
            std::cerr << "AFTER ALIGNMENT\n"
                      << "length of align row 0\t" << length(profileGaps) << "\n"
                      << "length of align row 1\t" << length(readGaps) << "\n"
                      << "\n"
                      << "ALIGNMENT\n"
                      << _profileGapsStructureStr(begin(profileGaps, Standard()), end(profileGaps, Standard())) << "\n"
                      << readGaps << "\n";

        if (options.debug)
        {
            swap(store.alignedReadStore, contigAlignedReads);
            std::cerr << "BEFORE INTEGRATION of " << store.readSeqStore[it->readId] << "\n";
            seqan::AlignedReadLayout layout;
            layoutAlignment(layout, store);
            std::stringstream ss;
            printAlignment(std::cerr, layout, store, front(store.alignedReadStore).contigId,
                           0, 1000, 0, 100);
            swap(store.alignedReadStore, contigAlignedReads);

            std::cerr << "begin/end positions of alignments after update\n";
            for (TAlignedReadIter it2 = begin(contigAlignedReads, Standard()); it2 != itEnd; ++it2)
            {
                typedef Gaps<TReadSeq, AnchorGaps<String<typename TFragmentStore::TReadGapAnchor> > > TReadGaps;
                TReadGaps readGaps(store.readSeqStore[it2->readId], it2->gaps);
                std::cerr << "it2->beginPos == " << it2->beginPos << ", it2->endPos == " << it2->endPos << "\n";
                std::cerr << "    " << source(readGaps) << "\n";
                if (windowBegin == windowEnd)
                    SEQAN_ASSERT_EQ((int)(it2->endPos - it2->beginPos), (int)length(readGaps));
            }
        }

        // Integrate back the profile and the new alignment.
        double beginIntegration = sysTime();
        if (options.debug)
            std::cerr << "WINDOW BEFORE UPDATE ALIGNMENTS [" << windowBegin << ", " << windowEnd << ")\n";
        _updateAlignments(windowBegin, windowEnd, profilePart, profileGaps, readGaps, info, windowInfo, it);
        SEQAN_ASSERT_LEQ(windowBegin, windowEnd);
        if (options.debug)
            std::cerr << "WINDOW AFTER UPDATE ALIGNMENTS [" << windowBegin << ", " << windowEnd << ")\n";
        times.realignIntegration = sysTime() - beginIntegration;

        if (options.debug)
        {
            std::cerr << "PROFILE PART AFTER ALIGNMENT UPDATE\n";
            _printProfile(std::cerr, profilePart);
        }

        replace(contigProfile, info.profBeginPos, info.profEndPos, profilePart);

        if (options.debug)
        {
            std::cerr << "PROFILE AFTER INTEGRATION PROFILE\n";
            _printProfile(std::cerr, contigProfile);

            swap(store.alignedReadStore, contigAlignedReads);
            std::cerr << "AFTER REALIGNMENT STEP of " << store.readSeqStore[it->readId] << "\n";
            seqan::AlignedReadLayout layout;
            layoutAlignment(layout, store);
            std::stringstream ss;
            printAlignment(std::cerr, layout, store, front(store.alignedReadStore).contigId,
                           0, 1000, 0, 100);
            swap(store.alignedReadStore, contigAlignedReads);

            std::cerr << "begin/end positions of alignments after update\n";
            for (TAlignedReadIter it2 = begin(contigAlignedReads, Standard()); it2 != itEnd; ++it2)
            {
                typedef Gaps<TReadSeq, AnchorGaps<String<typename TFragmentStore::TReadGapAnchor> > > TReadGaps;
                TReadGaps readGaps(store.readSeqStore[it2->readId], it2->gaps);
                std::cerr << "it2->beginPos == " << it2->beginPos << ", it2->endPos == " << it2->endPos << "\n";
                std::cerr << "    " << source(readGaps) << "\n";
                if (windowBegin == windowEnd)
                    SEQAN_ASSERT_EQ((int)(it2->endPos - it2->beginPos), (int)length(readGaps));
            }
        }
        SEQAN_ASSERT_LEQ(windowBegin, windowEnd);
    }

    // Register round's times.
    double t = sysTime() - roundStartTime;
    times.realignTotal += t;
    times.roundTimes.push_back(t);

    // Print times if configured to do so.
    if (options.printTiming)
    {
        std::cerr << "ANSON-MYERS REALIGNMENT TIMES (s)\n"
                  << "  main begin contig " << times.mainBeginContig << "\n"
                  << "  main end contig " << times.mainEndContig << "\n"
                  << "  main total " << times.mainTotal << "\n"
                  << "\n"
                  << "  realign extract profile " << times.realignExtractProfile << "\n"
                  << "  realign align profile " << times.realignAlign << "\n"
                  << "  realign integration profile " << times.realignIntegration << "\n"
                  << "  realign total " << times.realignTotal << "\n"
                  << "\n"
                  << "round times";
        for (unsigned i = 0; i < length(times.roundTimes); ++i)
            std::cerr << " " << times.roundTimes[i];
        std::cerr << "\n";
    }

    if (options.debug)
        std::cerr << "DONE ROUND WITH windowBegin == " << windowBegin << ", windowEnd == " << windowEnd << "\n";
}

// ----------------------------------------------------------------------------
// Function AnsonMyersRealignmentRound_::_updateAlignments()
// ----------------------------------------------------------------------------

// Return gaps/non-gaps string with c as the non-gaps character.
template <typename TGapsIter>
std::string _gapsStructureStr(TGapsIter begin, TGapsIter end)
{
    std::stringstream ss;
    for (; begin != end; ++begin)
        if (isGap(begin))
            ss << '-';
        else
            ss << *begin;
    return ss.str();
}

// Return gaps/non-gaps string with c as the non-gaps character.
template <typename TGapsIter>
std::string _profileGapsStructureStr(TGapsIter begin, TGapsIter end)
{
    // TGapsIter it = begin;
    std::stringstream ss;
    for (; begin != end; ++begin)
        if (isGap(begin))
        {
            ss << '-';
        }
        else
        {
            ProfileChar<Dna5> y = *begin;
            unsigned x = getMaxIndex(y);
            if (x == 5)
                ss << '=';
            else
                ss << Dna5(x);
        }

    // unsigned i = 0;
    // for (; it != end; ++it, ++i)
    // {
    //     for (unsigned x = 0; x < 5; ++x)
    //     {
    //         ProfileChar<Dna5> y = *it;
    //         // if (y.count[x])
    //         //     std::cout << "profile[" << i << "].count[ordValue(Dna5('" << Dna5(x) << "'))] = "
    //         //               << y.count[x] << ";\n";
    //         // if (y.count[5])
    //         //     std::cout << "profile[" << i << "].count[valueSize<Dna5>()] = " << y.count[5] << ";\n";
    //     }
    // }

    return ss.str();
}

template <typename TFragmentStore>
void AnsonMyersRealignmentRound_<TFragmentStore>::_updateAlignments(
        unsigned windowBegin,
        unsigned & windowEnd,  // only updated
        TProfileString & profilePart,
        // Align<TProfileString> /*const*/ & align,
        Gaps<TProfileString> /*const*/ & profileGaps,
        Gaps<TReadSeq> /*const*/ & readGaps,
        ExtractProfileInfo_ & info,
        WindowClippingInfo_ const & _windowInfo,
        TAlignedReadIter const & it)
{
    WindowClippingInfo_ windowInfo = _windowInfo;

    // typedef Align<TProfileString> TAlign;
    // typedef typename Row<TAlign>::Type TGaps;
    // typedef typename Iterator<TGaps, Standard>::Type TGapsIter;
    typedef Gaps<TProfileString> TProfileGaps;
    typedef typename Iterator<TProfileGaps>::Type TProfileGapsIter;
    typedef Gaps<TReadSeq> TReadGaps;
    typedef typename Iterator<TReadGaps>::Type TReadGapsIter;

    // Get copy of aligned read element.
    TAlignedReadElement el = *it, el2 = *it;

    // We build a new profile part.
    TProfileString newProfilePart;
    // reserve(newProfilePart, length(source(row(align, 0))));
    reserve(newProfilePart, length(source(profileGaps)));

    // Get iterators of alignment for profile (horizontal) and read (vertical).
    TProfileGapsIter itP = begin(profileGaps, Standard());
    TProfileGapsIter itPEnd = end(profileGaps, Standard());
    TReadGapsIter itR = begin(readGaps, Standard());
    TReadGapsIter itREnd = end(readGaps, Standard());
    // TGapsIter itP = begin(row(align, 0), Standard());
    // TGapsIter itPEnd = end(row(align, 0), Standard());
    // TGapsIter itR = begin(row(align, 1), Standard());
    // TGapsIter itREnd = end(row(align, 1), Standard());

    if (options.debug)
    {
        std::string str = _profileGapsStructureStr(itP, itPEnd);
        std::cerr << "PROFILE GAPS\t" << str << "\n"
                  << "READ GAPS   \t" << _gapsStructureStr(itR, itREnd) << "\n";
    }

    // Flags whether we are (possibly) in the initial or the trailing gap in the read/profile.
    bool inLeadingGapR = true, inLeadingGapP = true;
    if (windowInfo.clipLeft)
        inLeadingGapR = false;  // at least one base precedes current
    bool inTrailingGapR = false, inTrailingGapP = false;
    // Number of remaining alignment columns.
    // unsigned remainingCols = length(row(align, 0));
    unsigned remainingCols = length(profileGaps);

    // Reset begin/end position of alignment in profile
    info.aliBeginPos = info.aliEndPos = 0;

    // Construct read gaps, clear gaps from current alignment clipping (we get this through windowInfo) and keep a
    // counter for the position in the read gaps.  This is done by using a second copy since iterating and removing is
    // broken somehow.
    //
    // TODO(holtgrew): Iteration really broken after removal?
    typedef Gaps<TReadSeq, AnchorGaps<String<typename TFragmentStore::TReadGapAnchor> > > TAnchorReadGaps;
    TAnchorReadGaps anchorReadGaps(store.readSeqStore[el.readId], el.gaps);
    {
        TAnchorReadGaps anchorReadGaps2(store.readSeqStore[el2.readId], el2.gaps);
        if (options.debug)
            std::cerr << "removing gaps from " << windowInfo.readAliBeginPos << " to " << windowInfo.readAliEndPos
                      << "\n"
                      << "    from " << anchorReadGaps << "\n";
        clearGaps(anchorReadGaps);
        setClippedBeginPosition(anchorReadGaps, windowInfo.readEndPos);
        setClippedBeginPosition(anchorReadGaps2, windowInfo.readAliEndPos);
        if (options.debug)
            std::cerr << "copying gaps from >>" << anchorReadGaps2 << "<< to >>" << anchorReadGaps << "<<\n";
        copyGaps(anchorReadGaps, anchorReadGaps2);
        clearClipping(anchorReadGaps);
        clearClipping(anchorReadGaps2);

        setClippedEndPosition(anchorReadGaps, windowInfo.readBeginPos);
        setClippedEndPosition(anchorReadGaps2, windowInfo.readAliBeginPos);
        if (options.debug)
            std::cerr << "copying gaps from >>" << anchorReadGaps2 << "<< to >>" << anchorReadGaps << "<<\n";
        copyGaps(anchorReadGaps, anchorReadGaps2);
        clearClipping(anchorReadGaps);
        clearClipping(anchorReadGaps2);
        windowInfo.readAliEndPos -= (length(anchorReadGaps2) - length(anchorReadGaps));
        if (options.debug)
            std::cerr << "resulted in window on read " << windowInfo.readAliBeginPos << " to "
                      << windowInfo.readAliEndPos << "\n"
                      << "    for  " << anchorReadGaps << "\n";
    }
    int readGapsPos = windowInfo.readAliBeginPos;

    // Position in contig's multi-read alignment.
    unsigned pos = info.profBeginPos;

    // Iterate over both profile and read gaps and add the read bases back to the profile.  In case of new gaps, we have
    // to add a gap to the read alignments and also to the profile (as a new profile char).
    //
    // We record the positions of all-gaps columns and remove them later.  We also insert all-gaps columns into *it and
    // remove it together with the rest of the alignment.
    seqan::String<unsigned> gapPositions;
    for (; itP != itPEnd;)
    {
        SEQAN_ASSERT_NOT_MSG(isGap(itP) && isGap(itR), "No all-gaps columns in pairwise alignment.");

        // Without a countTrailingGaps() function, this is the best way to check whether we are in the trailing gap.
        inTrailingGapR = (inTrailingGapR || countGaps(itR) >= remainingCols);
        if (windowInfo.clipRight)
            inTrailingGapR = false;  // at least one base follows
        inTrailingGapP = (inTrailingGapP || countGaps(itP) >= remainingCols);

        if (isGap(itP))
        {
            // Case 1: gap in profile.
            //
            // Determine number of gaps in the profile/characters in the read and advance that far in profile and
            // multi-read alignment.  Extend the profile, adding the characters from the read alignment.  We need to
            // insert a gap in the whole multi-read alignment.
            unsigned numChars = std::min(countGaps(itP), countCharacters(itR));
            for (unsigned i = 0; i < numChars; ++i, ++itR, ++pos)
            {
                appendValue(newProfilePart, TProfileChar());  // new empty profile character
                TAlphabet x = *itR;
                back(newProfilePart).count[ordValue(x)] += 1;
                int numGaps = _insertGap(contigAlignedReads, pos, it);
                if (options.debug)
                    std::cerr << "inTrailingGapP == " << inTrailingGapP << "\n"
                              << "insertGap(contigAlignedRead, " << pos << ") == " << numGaps << "\n";
                if (windowBegin != windowEnd && !inTrailingGapP)
                    windowEnd += 1;  // gap inserted into profile/msa, increase window size
                back(newProfilePart).count[valueSize<TAlphabet>()] += numGaps;
                ++readGapsPos;
            }
            // Batch-update the profile gaps iterator and remaining column count.
            itP += numChars;
            remainingCols -= numChars;
            // We cannot be in the leading profile gap any more.
            inLeadingGapR = false;
            // Extend read alignment position.
            info.aliEndPos += numChars;
            if (options.debug)
                std::cerr << "info.aliEndPos += " << numChars << " == " << info.aliEndPos << "\n";
        }
        else if (isGap(itR))
        {
            // Case 2: gap in read (but not in profile, see assertion above).
            //
            // Determine number of gaps in the read/characters in the profile and advance that far in profile and
            // multi-read alignment.  Extend the profile, adding the gap from the read if we are neither in the leading
            // nore in the trailing part of the read's alignment row.
            unsigned numChars = std::min(countCharacters(itP), countGaps(itR));
            for (unsigned i = 0; i < numChars; ++i, ++itP, ++pos)
            {
                TProfileChar p = *itP;  // TODO(holtgrew): Copy because of proxy issues.
                appendValue(newProfilePart, p);
                if (!inLeadingGapR && !inTrailingGapR)
                {
                    back(newProfilePart).count[valueSize<TAlphabet>()] += 1;
                    if (options.debug)
                        std::cerr << "insertGap(anchorReadGaps, " << readGapsPos << ")\n";
                    insertGap(anchorReadGaps, readGapsPos++);
                    if (!inTrailingGapR)
                    {
                        info.aliEndPos += 1;
                        if (options.debug)
                            std::cerr << "info.aliEndPos += 1 == " << info.aliEndPos << "\n";
                    }
                }
                if (empty(p))
                    appendValue(gapPositions, pos);
            }
            // Batch-update the read gaps iterator and remaining column count.
            itR += numChars;
            remainingCols -= numChars;
            // We cannot be in the leading gap of the profile any more.
            inLeadingGapP = false;
            // If we are in the leading gap of the read, it begins (and ends) yet some character further to the right.
            if (inLeadingGapR)
            {
                info.aliBeginPos += numChars;
                info.aliEndPos += numChars;
                if (options.debug)
                    std::cerr << "info.aliBeginPos += " << numChars << " == " << info.aliBeginPos << "\n"
                              << "info.aliEndPos += " << numChars << " == " << info.aliEndPos << "\n";
            }
        }
        else
        {
            // Case 3: characters in both read and profile.
            //
            // Append the new character to the profile we create and update the counts.
            TProfileChar p = *itP;  // TODO(holtgrew): copy because of proxy issues
            appendValue(newProfilePart, p);
            TAlphabet x = *itR;
            back(newProfilePart).count[ordValue(x)] += 1;
            // Advance in profile and read row of alignment.
            ++itP;
            ++itR;
            // Advance in the read gaps and the position in the multi-read alignment.
            ++readGapsPos;
            ++pos;
            // One more character in read alignment and one fewer column to go.
            ++info.aliEndPos;
            --remainingCols;
            // Cannot be in leadin gaps any more.
            inLeadingGapR = inLeadingGapP = false;
        }
    }

    // Must be both at the end of the read and the profile alignment.
    SEQAN_ASSERT(itR == itREnd);
    SEQAN_ASSERT(itP == itPEnd);
    if (options.debug)
        std::cerr << "windowEnd == " << windowEnd << "\n";

    // -----------------------------------------------------------------------
    // Update the read alignment, write out new profile.
    // -----------------------------------------------------------------------

    // The read alignment in *it might have been updated in the process.  We update the position from info and the gaps
    // from el and info and windowInfo.

    if (options.debug)
        std::cerr << "REPLACING\tit->beginPos == " << it->beginPos << ", it->endPos == " << it->endPos << "\n";
    it->beginPos = info.profBeginPos + info.aliBeginPos - windowInfo.readAliBeginPos;
    it->endPos = it->beginPos + length(anchorReadGaps);
    swap(it->gaps, el.gaps);

    if (options.debug)
        std::cerr << "it->profBeginPos == " << info.profBeginPos << ", info.aliBeginPos == "
                  << info.aliBeginPos << "\n"
                  << "it->profEndPos == " << info.profEndPos << ", info.aliEndPos == " << info.aliEndPos << "\n"
                  << "windowInfo.readAliBeginPos == " << windowInfo.readAliBeginPos
                  << ", windowInfo.readAliEndPos == " << windowInfo.readAliEndPos << "\n"
                  << "it->beginPos == " << it->beginPos << ", it->endPos == " << it->endPos << "\n";
    int newContigProfileLength = length(contigProfile) + length(newProfilePart) - length(profilePart);
    (void)newContigProfileLength;  // only used for assertion
    SEQAN_ASSERT_LEQ(it->beginPos, newContigProfileLength);
    SEQAN_ASSERT_LEQ(it->endPos, newContigProfileLength);
    // SEQAN_ASSERT_EQ(it->endPos - it->beginPos, length(row(align, 0))); // TODO(holtgrew): Remove!

    // Write out the new profile.
    swap(newProfilePart, profilePart);

    // -----------------------------------------------------------------------
    // Remove all-gaps columns.
    // -----------------------------------------------------------------------

    if (options.debug && !empty(gapPositions))
        std::cerr << "removing all-gaps columns\n";
    for (unsigned i = 0; i < length(gapPositions); ++i)
    {
        if (options.debug)
            std::cerr << "removeGaps(contigAlignedReads, " << (gapPositions[i] - i)
                      << ") (info.profBeginPos == " << info.profBeginPos << ")\n";
        removeGap(contigAlignedReads, gapPositions[i] - i);
        SEQAN_ASSERT(empty(profilePart[gapPositions[i] - i - info.profBeginPos]));
        erase(profilePart, gapPositions[i] - i - info.profBeginPos);
        if (windowBegin != windowEnd)
            windowEnd -= 1;  // gap removed from profile/MSA, decrease window size
    }

#if SEQAN_ENABLE_TESTING
    {
        TAnchorReadGaps anchorReadGaps2(store.readSeqStore[it->readId], it->gaps);
        if (options.debug)
            std::cerr << "anchorReadGaps == " << anchorReadGaps2 << "\n";
        SEQAN_ASSERT_EQ(it->endPos - it->beginPos, (int)length(anchorReadGaps2));
    }
#endif  // #if SEQAN_ENABLE_TESTING
}

// ----------------------------------------------------------------------------
// Function reAlignment()
// ----------------------------------------------------------------------------

/*!
 * @fn reAlignment
 * @headerfile <seqan/realign.h>
 * @brief Perform realignment on a @link FragmentStore @endlink object.
 *
 * @signature void reAlignment(store, contigID, realignmentMethod, bandwidth, includeReference[,
 *                             windowBegin, windowEnd][, debug][, printTiming]);
 *
 * @param[in,out] store             The @link FragmentStore @endlink to perform realignment for.
 * @param[in]     realignmentMethod The realignment algorithm to use, <tt>1</tt> for affine gap costs, <tt>0</tt>
 *                                  for linear gap costs, affine gap costs are recommended, <tt>unsigned</tt>.
 * @param[in]     contigID          Identifier of the contig to realign, <tt>unsigned</tt>
 * @param[in]     bandwidth         Bandwidth to use in the pairwise DP alignment algorithms <tt>unsigned</tt>.
 * @param[in]     includeReference  A <tt>bool</tt>, if <tt>true</tt> then the reference will be included as a
 *                                  pseudo-read to increase stability in case of indels.  See section "Including
 *                                  the Reference" for details.
 * @param[in]     windowBegin       Optional <tt>unsigned</tt> window begin for windowed alignment,
 *                                  default: <tt>0</tt>, i.e. no windowing.  Also see section "Windowed
 *                                  Realignment" below.
 * @param[in]     windowEnd         Optional <tt>unsigned</tt> window end for windowed alignment, default: <tt>0</tt>,
 *                                  i.e. no windowing.  Also see section "Windowed Realignment" below.
 * @param[in]     debug             Optional <tt>bool</tt> to enable verbose logging, default: <tt>false</tt>.
 * @param[in]     printTiming       Optional <tt>bool</tt> to enable printing of times, default: <tt>false</tt>.
 *
 * This function implements a variant of the Anson-Myers algorithm to refine the multi-read alignments.  The algorithm
 * works in a round-robin fashion:
 *
 * First, the profile sequence is created for the multi-read alignment.  The algorithm then works in round.
 *
 * In each round, each read alignment is selected, removed from the profile and then aligned against the profile using
 * a DP algorithm.  For this alignment, a region of <tt>bandwidth / 2</tt> positions around the previous alignment
 * locus is chosen and aligned using a bandwidth of <tt>bandwidth</tt> in the DP alignment algorithm.  The profile is
 * then updated with the resulting pairwise alignment.  This is iterated until the global alignment score does not
 * improve.
 *
 * In the end, the reference sequence of the given contig is replaced by the consensus of the multi-read alignment.
 *
 * @section Including the Reference
 *
 * You can choose to include the reference sequence of the selected contig as a pseudo-read.  In this case, the
 * reference sequence of the selected contig will be added as an additional read and also a read alignment will be
 * created for it.  The appended pseudo-read will have the highest read ID in the store.
 *
 * @section Windowed Realignment
 *
 * Sometimes, it is useful to just realign a part of the fragment store.  This can be done by using the parameters
 * <tt>windowBegin</tt> and <tt>windowEnd</tt>.  This refers to view positions in the current multi-read alignment.
 * Set both to <tt>0</tt> (also the default) to disalbe windowed realignment.
 *
 * @section References
 *
 * <ul><li>Anson, Eric L., and Eugene W. Myers. "ReAligner: a program for refining DNA sequence multi-alignments."
 *         Journal of Computational Biology 4.3 (1997): 369-383.</li></ul>
 */

template <typename TSpec, typename TConfig>
void reAlignment(FragmentStore<TSpec, TConfig> & store,
                 unsigned contigID,
                 unsigned realignmentMethod,
                 unsigned bandwidth,
                 bool includeReference,
                 unsigned windowBegin = 0,
                 unsigned windowEnd = 0,
                 bool debug = false,
                 bool printTiming = false)
{
    RealignmentOptions_ options;
    options.method = realignmentMethod ? RealignmentOptions_::ANSON_MYERS_GOTOH : RealignmentOptions_::ANSON_MYERS_NW;
    options.bandwidth = bandwidth;
    options.environment = bandwidth / 2;
    options.includeReference = includeReference;
    options.debug = debug;
    options.printTiming = printTiming;

    AnsonMyersRealigner_<FragmentStore<TSpec, TConfig> > realigner(store, options);
    realigner.run(contigID, windowBegin, windowEnd);
}

template <typename TSpec, typename TConfig, typename TScore>
void reAlignment(FragmentStore<TSpec, TConfig> & store,
                 TScore & score,
                 unsigned contigID,
                 unsigned realignmentMethod,
                 unsigned bandwidth,
                 bool includeReference,
                 unsigned windowBegin = 0,
                 unsigned windowEnd = 0,
                 bool debug = false,
                 bool printTiming = false)
{
    (void)score;
    reAlignment(store, contigID, realignmentMethod, bandwidth, includeReference, windowBegin, windowEnd, debug,
                printTiming);
}

}  // namespace seqan

#endif  // INCLUDE_SEQAN_REALIGN_REALIGN_BASE_H_
