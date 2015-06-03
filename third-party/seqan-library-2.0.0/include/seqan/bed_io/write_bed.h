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
// Writing of BED to files.
// ==========================================================================

// TODO(holtgrew): Add overload with BedIOContext.

#ifndef INCLUDE_SEQAN_BED_IO_WRITE_BED_H_
#define INCLUDE_SEQAN_BED_IO_WRITE_BED_H_

namespace seqan {

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function writeRecord()                                      [Bed3 BedRecord]
// ----------------------------------------------------------------------------

// Similar to readRecord() for BedRecord, we have various overloads of _writeBedRecord().

template <typename TTarget>
inline void
_writeBedRecord(TTarget & target, BedRecord<Bed3> const & record, CharString const & ref)
{
    write(target, ref);
    writeValue(target, '\t');
    // NB: in contrast to many other text-based formats, UCSC BED uses 0-based and not 1-based coordinates.
    appendNumber(target, record.beginPos);
    writeValue(target, '\t');
    appendNumber(target, record.endPos);
}

template <typename TTarget>
inline void
_writeBedRecord(TTarget & target, BedRecord<Bed4> const & record, CharString const & ref)
{
    _writeBedRecord(target, static_cast<BedRecord<Bed3> const &>(record), ref);
    writeValue(target, '\t');
    write(target, record.name);
}

template <typename TTarget>
inline void
_writeBedRecord(TTarget & target, BedRecord<Bed5> const & record, CharString const & ref)
{
    _writeBedRecord(target, static_cast<BedRecord<Bed4> const &>(record), ref);
    writeValue(target, '\t');
    write(target, record.score);
}

template <typename TTarget>
inline void
_writeBedRecord(TTarget & target, BedRecord<Bed6> const & record, CharString const & ref)
{
    _writeBedRecord(target, static_cast<BedRecord<Bed5> const &>(record), ref);
    writeValue(target, '\t');
    writeValue(target, record.strand);
}

template <typename TTarget>
inline void
_writeBedRecord(TTarget & target, BedRecord<Bed12> const & record, CharString const & ref)
{
    _writeBedRecord(target, static_cast<BedRecord<Bed6> const &>(record), ref);
    writeValue(target, '\t');

    appendNumber(target, record.thickBegin);
    writeValue(target, '\t');

    appendNumber(target, record.thickEnd);
    writeValue(target, '\t');

    appendNumber(target, record.itemRgb.red);
    writeValue(target, ',');
    appendNumber(target, record.itemRgb.green);
    writeValue(target, ',');
    appendNumber(target, record.itemRgb.blue);
    writeValue(target, '\t');

    appendNumber(target, record.blockCount);
    writeValue(target, '\t');

    for (unsigned i = 0; i < length(record.blockSizes); ++i)
    {
        if (i > 0)
            writeValue(target, ',');
        appendNumber(target, record.blockSizes[i]);
    }

    writeValue(target, '\t');

    for (unsigned i = 0; i < length(record.blockBegins); ++i)
    {
        if (i > 0)
            writeValue(target, ',');
        appendNumber(target, record.blockBegins[i]);
    }
}

template <typename TTarget, typename TRecordSpec>
inline void
writeRecord(TTarget & target, BedRecord<TRecordSpec> const & record, Bed const & /*tag*/)
{
    _writeBedRecord(target, record, record.ref);
    if (!empty(record.data))
    {
        writeValue(target, '\t');
        write(target, record.data);
    }
    writeValue(target, '\n');
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_BED_IO_WRITE_BED_H_
