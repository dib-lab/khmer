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

#ifndef INCLUDE_SEQAN_ROI_IO_ROI_RECORD_H_
#define INCLUDE_SEQAN_ROI_IO_ROI_RECORD_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class RoiHeader
// ----------------------------------------------------------------------------

class RoiHeader
{
public:
    seqan::StringSet<seqan::CharString> extraColumns;
};

// ----------------------------------------------------------------------------
// Class RoiRecord
// ----------------------------------------------------------------------------

// One entry in a ROI file.

class RoiRecord
{
public:
    static const int INVALID_POS = -1;

    // Chromosome that this region is on.
    CharString ref;

    // Begin and end position.
    int beginPos;
    int endPos;

    // The strand of this region.
    char strand;  // TODO(holtgrew): Rename to strand and make char?

    // The name of the region.
    CharString name;

    // The coverage over the length of the region.
    String<unsigned> count;  // TODO(holtgrew): Rename to coverages/counts?

    // The length of the region.
    unsigned len;

    // The largest count.
    unsigned countMax;  // TODO(holtgrew): Rename to maxCount/maxCoverages?

    // Additional data as string.
    seqan::StringSet<seqan::CharString> data;

    RoiRecord() : beginPos(INVALID_POS), endPos(INVALID_POS), strand('.'), len(0), countMax(0)
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

inline void clear(RoiRecord & record)
{
    clear(record.ref);
    clear(record.count);
    clear(record.data);
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ROI_IO_ROI_RECORD_H_
