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

#ifndef INCLUDE_SEQAN_REALIGN_REALIGN_HELPER_FUNCTIONS_H_
#define INCLUDE_SEQAN_REALIGN_REALIGN_HELPER_FUNCTIONS_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TPos, typename TGapAnchor, typename TSpec, typename TGapPos>
inline bool  // true if removed gap that is not leading/trailing gap
_removeGap2(AlignedReadStoreElement<TPos, TGapAnchor, TSpec> & alignedRead,
            TGapPos const gapPos)
{
    typedef String<TGapAnchor> TGaps;
    typedef typename Iterator<TGaps, Standard>::Type TGapIter;
    if (gapPos < (TGapPos) alignedRead.beginPos)
    {
        --alignedRead.beginPos; --alignedRead.endPos;
    }
    else if (gapPos < (TGapPos) alignedRead.endPos)
    {
        --alignedRead.endPos;
        TGapIter gapIt = upperBoundGapAnchor(alignedRead.gaps, gapPos - alignedRead.beginPos, SortGapPos());
        TGapIter gapItEnd = end(alignedRead.gaps, Standard());
        // Note: We might create empty gaps here
        for(; gapIt != gapItEnd; ++gapIt)
            --(gapIt->gapPos);
        return true;
    }
    return false;
}

// Returns number of returned gaps.
template <typename TAlignedReads, typename TSpec, typename TGapPos, typename TAlignedReadIt>
inline int
_removeGap(String<TAlignedReads, TSpec>& alignedReadStore,
           TGapPos const gapPos,
           TAlignedReadIt const & skipIt)
{
    typedef String<TAlignedReads, TSpec> TAlignedReadStore;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type TAlignIter;
    int result = 0;

    TAlignIter alignIt = begin(alignedReadStore, Standard());
    TAlignIter alignItEnd = end(alignedReadStore, Standard());
    for(;alignIt != alignItEnd; ++alignIt)
        if (alignIt != skipIt)
            result += _removeGap2(*alignIt, gapPos);
    return result;
}

template <typename TAlignedReads, typename TSpec, typename TGapPos, typename TAlignedReadIt>
inline int
_insertGap(String<TAlignedReads, TSpec>& alignedReadStore,
           TGapPos const gapPos,
           TAlignedReadIt const & skipIt)
{
    typedef String<TAlignedReads, TSpec> TAlignedReadStore;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type TAlignIter;

    int numGaps = 0;
    TAlignIter alignIt = begin(alignedReadStore, Standard());
    TAlignIter alignItEnd = end(alignedReadStore, Standard());
    for(;alignIt != alignItEnd; ++alignIt)
        if (alignIt != skipIt)
            numGaps += insertGap(*alignIt, gapPos);
    return numGaps;
}

// ----------------------------------------------------------------------------
// Helper Function _printProfile()
// ----------------------------------------------------------------------------

// 0   - space
// 1-9 - 1..9
// >10 - '!' + ceil(log(count))
inline char _countToLogChar(int count)
{
    if (count <= 0)
        return ' ';
    else if (count <= 9)
        return '0' + count;
    else
        return '!' + int(ceil(log10(1.0 * count)));
}

template <typename TStream, typename TAlphabet, typename TSpec>
void _printProfile(TStream & stream,
                   String<ProfileChar<TAlphabet>, TSpec> const & profile)
{
    stream << "PROFILE (0 - space, 1..9 - 1..9, >10 - ! + int(ceil(log10(count)))\n"
           << "------------------------------------------------------------------\n";

    stream << "   ";
    for (unsigned i = 0 ; i < length(profile); ++i)
    {
        if (i % 10 == 0)
            stream << (i / 100) % 10;
        else if (i % 5 == 0)
            stream << (i / 100) % 10;
        else
            stream << ' ';
    }
    stream << "\n";

    stream << "   ";
    for (unsigned i = 0 ; i < length(profile); ++i)
    {
        if (i % 10 == 0)
            stream << (i / 10) % 10;
        else if (i % 5 == 0)
            stream << (i / 10) % 10;
        else
            stream << ' ';
    }
    stream << "\n";

    stream << "   ";
    for (unsigned i = 0 ; i < length(profile); ++i)
    {
        if (i % 10 == 0)
            stream << '0';
        else if (i % 5 == 0)
            stream << '5';
        else
            stream << ' ';
    }
    stream << "\n";

    stream << "   ";
    for (unsigned i = 0 ; i < length(profile); ++i)
    {
        if (i % 10 == 0)
            stream << ':';
        else if (i % 5 == 0)
            stream << '.';
        else
            stream << ' ';
    }
    stream << "\n";

    for (unsigned i = 0; i <= valueSize<TAlphabet>(); ++i)
    {
        if (i == valueSize<TAlphabet>())
            stream << "-: ";
        else
            stream << TAlphabet(i) << ": ";

        for (unsigned j = 0; j < length(profile); ++j)
            stream << _countToLogChar(profile[j].count[i]);

        stream << "\n";
    }
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_REALIGN_REALIGN_HELPER_FUNCTIONS_H_
