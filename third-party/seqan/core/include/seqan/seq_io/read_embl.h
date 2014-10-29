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
// Read support for the EMBL file format.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_SEQ_IO_READ_EMBL_H_
#define EXTRAS_INCLUDE_SEQAN_SEQ_IO_READ_EMBL_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*
struct Embl_;
typedef Tag<Embl_> Embl;
*/

struct EmblHeader_;
typedef Tag<EmblHeader_> EmblHeader;

struct EmblSequence_;
typedef Tag<EmblSequence_> EmblSequence;

enum EmblErrorCodes_
{
    IOERR_EMBL_WRONG_RECORD = 2048
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function splitEmblHeader()
// ----------------------------------------------------------------------------

/**
.Function.splitEmblHeader
..cat:Input/Output
..signature:startsWith(key, value, line)
..summary:Split an EMBL header line.
..param.key:The 2-character header type.
..param.value:The line's value.
..param.line:The header line to split.
..returns:$void$
..include:seqan/stream.h
*/

template <typename TKey, typename TValue, typename TLine>
inline void
splitEmblHeader(TKey & key, TValue & value, TLine const & line)
{
    clear(key);

    enum State { IN_KEY, IN_SPACE, IN_VALUE };
    State state = IN_KEY;

    typedef typename Iterator<TLine const, Rooted>::Type TIterator;
    TIterator it = begin(line, Rooted());
    for (; !atEnd(it); goNext(it))
    {
        if (state == IN_KEY)
        {
            appendValue(key, *it);
            if (length(key) == 2u)
                state = IN_SPACE;
        }
        else if (state == IN_SPACE && !isblank(*it))
        {
            break;
        }
    }

    value = suffix(line, position(it));
}

// ----------------------------------------------------------------------------
// Function nextIs()
// ----------------------------------------------------------------------------

// nextIs() for EMBL header records.

template <typename TStream, typename TSpec>
inline bool
nextIs(RecordReader<TStream, TSpec> & reader, EmblHeader const & /*tag*/)
{
    return !atEnd(reader) && !isblank(value(reader));
}

// nextIs() for EMBL sequence records.

template <typename TStream, typename TSpec>
inline bool
nextIs(RecordReader<TStream, TSpec> & reader, EmblSequence const & /*tag*/)
{
    return !atEnd(reader) && isblank(value(reader));
}

// ----------------------------------------------------------------------------
// Function readRecord()
// ----------------------------------------------------------------------------

// readRecord() for EMBL header records.

template <typename TTarget, typename TStream, typename TSpec>
inline int
readRecord(TTarget & buffer, RecordReader<TStream, SinglePass<TSpec> > & reader,
           EmblHeader const & tag)
{
    if (atEnd(reader))
        return EOF_BEFORE_SUCCESS;
    if (!nextIs(reader, tag))
        return IOERR_EMBL_WRONG_RECORD;
    clear(buffer);

    return readLine(buffer, reader);
}

// readRecord() for EMBL sequences records.

// Read all sequence, eat/ignore '//' line.

template <typename TTarget, typename TStream, typename TSpec>
inline int
readRecord(TTarget & buffer, RecordReader<TStream, SinglePass<TSpec> > & reader,
           EmblSequence const & tag)
{
    if (atEnd(reader))
        return EOF_BEFORE_SUCCESS;
    if (!nextIs(reader, tag))
        return IOERR_EMBL_WRONG_RECORD;
    clear(buffer);

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
        appendValue(buffer, value(reader));
    }

    return EOF_BEFORE_SUCCESS;
}

// readRecord() for EMBL id/seq pairs.

template <typename TId, typename TSequence, typename TStream, typename TSpec>
inline int
readRecord(TId & id,
           TSequence & sequence,
           RecordReader<TStream, SinglePass<TSpec> > & reader,
           Embl const & /*tag*/)
{
    CharString buffer;
    while (nextIs(reader, EmblHeader()))
    {
        int res = readRecord(buffer, reader, EmblHeader());
        if (res != 0)
            return res;
        if (startsWith(buffer, "ID"))
        {
            CharString k, v;
            splitEmblHeader(k, id, buffer);
        }
    }
    clear(sequence);
    int res = readRecord(sequence, reader, EmblSequence());
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
          Embl const & tag)
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

#endif  // #ifndef EXTRAS_INCLUDE_SEQAN_SEQ_IO_READ_EMBL_H_
