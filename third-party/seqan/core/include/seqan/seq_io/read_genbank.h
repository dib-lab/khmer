// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Read support for the GenBank file format.
//
// Note: We do not want to require a lookahead because then we might need too
// many buffers.  Thus, we heuristically only look at the first character of
// the field name to differentiate between headers and the sequence start
// label "ORIGIN".
//
// The sequence identifier for record-reading is the VERSION field.
// ==========================================================================


#ifndef EXTRAS_INCLUDE_SEQAN_SEQ_IO_READ_GENBANK_H_
#define EXTRAS_INCLUDE_SEQAN_SEQ_IO_READ_GENBANK_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct GenBank_;
typedef Tag<GenBank_> GenBank;

struct GenBankHeader_;
typedef Tag<GenBankHeader_> GenBankHeader;

struct GenBankSequence_;
typedef Tag<GenBankSequence_> GenBankSequence;

enum GenBankErrorCodes_
{
    IOERR_GENBANK_WRONG_RECORD = 2048
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function splitGenBankHeader()
// ----------------------------------------------------------------------------

/**
.Function.splitGenBankHeader
..cat:Input/Output
..signature:splitGenBankHeader(key, value, lines)
..summary:Split an GenBank header field/value.
..remarks:You can only call this on a whole top-level field such as $SOURCE$, possibly including the following subfields! You cannot split out subfield values.
..param.key:The 2-character header type.
..param.value:The line's value.
..param.line:The lines with the header field to split.
..returns:$void$
..include:seqan/stream.h
..see:Function.splitEmblHeader
*/

template <typename TKey, typename TValue, typename TLine>
inline void
splitGenBankHeader(TKey & key, TValue & value, TLine const & lines)
{
    splitEmblHeader(key, value, lines);
    clear(key);

    enum State { IN_KEY, IN_SPACE, IN_VALUE };
    State state = IN_KEY;

    typedef typename Iterator<TLine const, Rooted>::Type TIterator;
    TIterator it = begin(lines, Rooted());
    for (; !atEnd(it); goNext(it))
    {
        if (state == IN_KEY)
        {
            if (isblank(*it))
                state = IN_SPACE;
            else
                appendValue(key, *it);
        }
        else if (state == IN_SPACE && !isblank(*it))
        {
            break;
        }
    }

    value = suffix(lines, position(it));
}

// ----------------------------------------------------------------------------
// Function nextIs()
// ----------------------------------------------------------------------------

// nextIs() for GenBank header records.

template <typename TStream, typename TSpec>
inline bool
nextIs(RecordReader<TStream, TSpec> & reader, GenBankHeader const & /*tag*/)
{
    return !atEnd(reader) && value(reader) != 'O' && !isblank(value(reader));
}

// nextIs() for GenBank sequence records.

template <typename TStream, typename TSpec>
inline bool
nextIs(RecordReader<TStream, TSpec> & reader, GenBankSequence const & /*tag*/)
{
    return !atEnd(reader) && value(reader) == 'O';
}

// ----------------------------------------------------------------------------
// Function readRecord()
// ----------------------------------------------------------------------------

// readRecord() for GenBank header records.

// Normalize in-field line endings to '\n'.

template <typename TTarget, typename TStream, typename TSpec>
inline int
readRecord(TTarget & result, RecordReader<TStream, SinglePass<TSpec> > & reader,
           GenBankHeader const & tag)
{
    if (atEnd(reader))
        return EOF_BEFORE_SUCCESS;
    if (!nextIs(reader, tag))
        return IOERR_GENBANK_WRONG_RECORD;
    clear(result);

    int res = readLine(result, reader);
    if (res != 0)
        return res;
    while (!atEnd(reader) && isblank(value(reader)))
    {
        appendValue(result, '\n');
        res = readLine(result, reader);
        if (res != 0)
            return res;
    }

    return 0;
}

// readRecord() for GenBank sequences records.

// Read all sequence, eat/ignore '//' line.

template <typename TTarget, typename TStream, typename TSpec>
inline int
readRecord(TTarget & result, RecordReader<TStream, SinglePass<TSpec> > & reader,
           GenBankSequence const & tag)
{
    if (atEnd(reader))
        return EOF_BEFORE_SUCCESS;
    if (!nextIs(reader, tag))
        return IOERR_GENBANK_WRONG_RECORD;
    clear(result);

    // Skip 'ORIGIN' line.
    CharString buffer;
    int res = readLine(buffer, reader);
    if (res != 0)
        return res;
    if (!startsWith(buffer, "ORIGIN"))
        return IOERR_GENBANK_WRONG_RECORD;

    // Read sequence.
    for (; !atEnd(reader); goNext(reader))
    {
        if (isspace(value(reader)) || isdigit(value(reader)))
            continue;  // Skip blanks.
        // If we found a slash, look for the next one.
        if (value(reader) == '/')
        {
            if (goNext(reader))  // Returns true if at end.
                break;
            if (value(reader) == '/')
            {
                skipLine(reader);
                return 0;  // OK, found complete terminator.
            }
        }

        // Otherwise, append the just found character to buffer.
        appendValue(result, value(reader));
    }

    return EOF_BEFORE_SUCCESS;
}

// readRecord() for GenBank id/seq pairs.

template <typename TId, typename TSequence, typename TStream, typename TSpec>
inline int
readRecord(TId & id,
           TSequence & sequence,
           RecordReader<TStream, SinglePass<TSpec> > & reader,
           GenBank const & /*tag*/)
{
    CharString buffer;
    while (nextIs(reader, GenBankHeader()))
    {
        int res = readRecord(buffer, reader, GenBankHeader());
        if (res != 0)
            return res;
        if (startsWith(buffer, "VERSION"))
        {
            CharString k;
            splitGenBankHeader(k, id, buffer);
        }
    }
    clear(sequence);
    int res = readRecord(sequence, reader, GenBankSequence());
    if (res != 0)
        return res;
    
    return 0;
}

// ----------------------------------------------------------------------------
// Function read2()
// ----------------------------------------------------------------------------

template <typename TIdString, typename TIdSpec, typename TSeqString, typename TSeqSpec,
          typename TStream, typename TSpec>
int read2(StringSet<TIdString, TIdSpec> & sequenceIds,
          StringSet<TSeqString, TSeqSpec> & sequences,
          RecordReader<TStream, SinglePass<TSpec> > & reader,
          GenBank const & tag)
{
    TIdString id;
    TSeqString seq;

    while (!atEnd(reader))
    {
        int res = readRecord(id, seq, reader, tag);
        if (res != 0)
            return res;
        appendValue(sequenceIds, id);
        appendValue(sequences, seq);
    }

    return 0;
}

}  // namespace seqan

#endif  // #ifndef EXTRAS_INCLUDE_SEQAN_SEQ_IO_READ_GENBANK_H_
