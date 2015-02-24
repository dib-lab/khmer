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

#ifndef INCLUDE_SEQAN_ROI_IO_READ_ROI_H_
#define INCLUDE_SEQAN_ROI_IO_READ_ROI_H_

#include <seqan/stream.h>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

struct Roi_;
typedef Tag<Roi_> Roi;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Class MagicHeader
// ----------------------------------------------------------------------------

template <typename T>
struct MagicHeader<Roi, T> :
    public MagicHeader<Nothing, T> {};

// ----------------------------------------------------------------------------
// Class FileExtensions
// ----------------------------------------------------------------------------

template <typename T>
struct FileExtensions<Roi, T>
{
    static char const * VALUE[1];    // default is one extension
};

template <typename T>
char const * FileExtensions<Roi, T>::VALUE[1] =
{
    ".roi"     // default output extension
};


// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord()                                            [RoiHeader]
// ----------------------------------------------------------------------------

template <typename TForwardIter>
void readHeader(RoiHeader & header, RoiIOContext &, TForwardIter & iter, Roi const & /*tag*/)
{
    typedef OrFunctor<IsTab, IsNewline> TNextEntry;

    // skip first hash
    skipOne(iter, EqualsChar<'#'>());
    if (!atEnd(iter) && *iter != '#')
        skipLine(iter);  // skip ROI line if any.

    // check if column header present
    if (!atEnd(iter) && *iter != '#')
        return;  // no header
    // skip run of hashes
    skipUntil(iter, NotFunctor<EqualsChar<'#'> >());

    // skip REF\t
    skipUntil(iter, TNextEntry());
    skipOne(iter, EqualsChar<'\t'>());

    // skip START\t
    skipUntil(iter, TNextEntry());
    skipOne(iter, EqualsChar<'\t'>());

    // skip END\t
    skipUntil(iter, TNextEntry());
    skipOne(iter, EqualsChar<'\t'>());

    // skip NAME\t
    skipUntil(iter, TNextEntry());
    skipOne(iter, EqualsChar<'\t'>());

    // skip LENGTH\t
    skipUntil(iter, TNextEntry());
    skipOne(iter, EqualsChar<'\t'>());

    // skip STRAND\t
    skipUntil(iter, TNextEntry());
    skipOne(iter, EqualsChar<'\t'>());

    // skip MAX_COUNT
    skipUntil(iter, TNextEntry());

    do
    {
        resize(header.extraColumns, length(header.extraColumns) + 1);
        skipOne(iter, EqualsChar<'\t'>());
        readUntil(header.extraColumns[length(header.extraColumns) - 1], iter, TNextEntry());
    }
    while (!atEnd(iter) && !IsNewline()(*iter));

    resize(header.extraColumns, length(header.extraColumns) - 1);
    skipLine(iter);
}

// ----------------------------------------------------------------------------
// Function readRecord()                                            [RoiRecord]
// ----------------------------------------------------------------------------

template <typename TForwardIter>
void readRecord(RoiRecord & record, RoiIOContext & context, TForwardIter & iter, Roi const & /*tag*/)
{
    typedef OrFunctor<IsTab, IsNewline> TNextEntry;

    // Read reference name.
    clear(record.ref);
    readUntil(record.ref, iter, TNextEntry());

    // Skip TAB.
    skipOne(iter, IsTab());

    // Read and parse start position.
    clear(context.buffer);
    readUntil(context.buffer, iter, TNextEntry());
    if (!lexicalCast(record.beginPos, context.buffer))
        throw BadLexicalCast(record.beginPos, context.buffer);
    record.beginPos -= 1;  // transform to 0-based

    // Skip TAB.
    skipOne(iter, IsTab());

    // Read and parse end position.
    clear(context.buffer);
    readUntil(context.buffer, iter, TNextEntry());
    if (!lexicalCast(record.endPos, context.buffer))
        throw BadLexicalCast(record.endPos, context.buffer);

    // Skip TAB.
    skipOne(iter, IsTab());

    // Read and parse region name.
    clear(record.name);
    readUntil(record.name, iter, TNextEntry());

    // Skip TAB.
    skipOne(iter, IsTab());

    // Read and parse region length.
    clear(context.buffer);
    readUntil(context.buffer, iter, TNextEntry());
    if (!lexicalCast(record.len, context.buffer))
        throw BadLexicalCast(record.len, context.buffer);

    // Skip TAB.
    skipOne(iter, IsTab());

    // Read strand.
    readOne(record.strand, iter, OrFunctor<OrFunctor<EqualsChar<'.'>, EqualsChar<'+'> >, EqualsChar<'-'> >());

    // Skip TAB.
    skipOne(iter, IsTab());

    // Read max count.
    clear(context.buffer);
    readUntil(context.buffer, iter, TNextEntry());
    if (!lexicalCast(record.countMax, context.buffer))
        throw BadLexicalCast(record.countMax, context.buffer);

    // Skip TAB.
    skipOne(iter, IsTab());

    // Read data field.
    do
    {
        clear(context.buffer);
        readUntil(context.buffer, iter, TNextEntry());
        if (!atEnd(iter))
        {
            if (!IsNewline()(value(iter)))
                appendValue(record.data, context.buffer);
            else
                break;

            skipOne(iter);
        }
    }
    while (true);
    skipLine(iter);

    // Individual counts.
    clear(record.count);
    DirectionIterator<String<char>, Input>::Type castIter = begin(context.buffer);

    while (!atEnd(castIter))
    {
        clear(context.castBuffer);
        readUntil(context.castBuffer, castIter, OrFunctor<EqualsChar<','>, IsNewline>());
        if (!empty(context.castBuffer))
        {
            unsigned count = 0;
            if (!lexicalCast(count, context.castBuffer))
                throw BadLexicalCast(count, context.castBuffer);
            appendValue(record.count, count);
            record.countMax = std::max(record.countMax, back(record.count));
        }
        if (!atEnd(castIter))
            skipOne(castIter);
    }
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ROI_IO_READ_ROI_H_
