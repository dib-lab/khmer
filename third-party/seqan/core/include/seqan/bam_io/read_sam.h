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
// Code for reading SAM.
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_READ_SAM_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_READ_SAM_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Tag.Sam
..cat:BAM I/O
..signature:Sam
..summary:Tag for identifying the SAM format.
..include:seqan/bam_io.h
..see:Tag.Bam
*/

struct Sam_;
typedef Tag<Sam_> const Sam;

enum SamTokenizeErrors_
{
    SAM_INVALID_RECORD = 2048
};

struct SamHeader_;
typedef Tag<SamHeader_> SamHeader;

struct SamAlignment_;
typedef Tag<SamAlignment_> SamAlignment;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function nextIs()                                                  SamHeader
// ----------------------------------------------------------------------------

template <typename TStream, typename TPass>
inline bool nextIs(RecordReader<TStream, TPass> & reader, SamHeader const & /*tag*/)
{
    if (atEnd(reader))
        return false;
    return value(reader) == '@';
}

// ----------------------------------------------------------------------------
// Function nextIs()                                               SamAlignment
// ----------------------------------------------------------------------------

template <typename TStream, typename TPass>
inline bool nextIs(RecordReader<TStream, TPass> & reader, SamAlignment const & /*tag*/)
{
    if (atEnd(reader))
        return false;
    return value(reader) != '@';
}

// ----------------------------------------------------------------------------
// Function skipRecord()                                              SamHeader
// ----------------------------------------------------------------------------

template <typename TStream, typename TPass>
inline int skipRecord(RecordReader<TStream, TPass> & reader,
                      SamHeader const & tag)
{
    if (atEnd(reader))
        return EOF_BEFORE_SUCCESS;
    if (!nextIs(reader, tag))
        return SAM_INVALID_RECORD;
    int res = skipLine(reader);
    if (res == 0 || res == EOF_BEFORE_SUCCESS)
        return 0;
    else
        return res;
}

// ----------------------------------------------------------------------------
// Function skipRecord()                                           SamAlignment
// ----------------------------------------------------------------------------

template <typename TStream, typename TPass>
inline int skipRecord(RecordReader<TStream, TPass> & reader,
                      SamAlignment const & tag)
{
    if (atEnd(reader))
        return EOF_BEFORE_SUCCESS;
    if (!nextIs(reader, tag))
        return SAM_INVALID_RECORD;
    int res = skipLine(reader);
    if (res == 0 || res == EOF_BEFORE_SUCCESS)
        return 0;
    else
        return res;
}

// ----------------------------------------------------------------------------
// Function readRecord()                                        BamHeaderRecord
// ----------------------------------------------------------------------------

template <typename TStream, typename TSpec, typename TNameStore, typename TNameStoreCache>
int readRecord(BamHeaderRecord & record,
               BamIOContext<TNameStore, TNameStoreCache> & context,
               RecordReader<TStream, SinglePass<TSpec> > & reader,
               Sam const & /*tag*/)
{
    clear(record);

    // Make sure the first character is '@'.
    char c = value(reader);
    if (c != '@')
        return SAM_INVALID_RECORD;
    if (goNext(reader))
        return SAM_INVALID_RECORD;

    // Read the header tag.
    char c1 = value(reader);
    if (goNext(reader))
        return SAM_INVALID_RECORD;
    char c2 = value(reader);
    if (goNext(reader))
        return SAM_INVALID_RECORD;
    // Determine header type.
    if (c1 == 'H' && c2 == 'D')
        record.type = BAM_HEADER_FIRST;
    else if (c1 == 'S' && c2 == 'Q')
        record.type = BAM_HEADER_REFERENCE;
    else if (c1 == 'R' && c2 == 'G')
        record.type = BAM_HEADER_READ_GROUP;
    else if (c1 == 'P' && c2 == 'G')
        record.type = BAM_HEADER_PROGRAM;
    else if (c1 == 'C' && c2 == 'O')
        record.type = BAM_HEADER_COMMENT;
    else
        return SAM_INVALID_RECORD;

    if (record.type == BAM_HEADER_COMMENT)
    {
        int res = skipChar(reader, '\t');
        if (res != 0)
            return res;
        CharString &buffer = context.buffer;
        res = readLine(buffer, reader);
        if (res != 0)
            return res;
        appendValue(record.tags, Pair<CharString>(CharString(), buffer));
    }
    else
    {
        // Read the rest of the line into the tag field of record.
        int res = 0;
        CharString key, val;
        while (!atEnd(reader) && value(reader) == '\t')
        {
            clear(key);
            clear(val);

            res = skipChar(reader, '\t');
            if (res != 0)
                return res;
            res = readUntilChar(key, reader, ':');
            if (res != 0)
                return res;
            if (goNext(reader))
                return SAM_INVALID_RECORD;
            res = readUntilOneOf(val, reader, '\t', '\r', '\n');
            if (res != 0 && res != EOF_BEFORE_SUCCESS)
                return res;

            appendValue(record.tags, Pair<CharString>(key, val));
        }
    }

    // Skip remaining line break.
    int res = skipLine(reader);
    if (res != 0 && res != EOF_BEFORE_SUCCESS)
        return res;
    return 0;
}

// ----------------------------------------------------------------------------
// Function readRecord()                                              BamHeader
// ----------------------------------------------------------------------------

/**
.Function.readRecord
..signature:readRecord(headerRecord, context, recordReader, tag)
..param.recordReader:The RecordReader to read from.
...type:Class.RecordReader
...remarks:Use for SAM.
*/

template <typename TStream, typename TSpec, typename TNameStore, typename TNameStoreCache>
int readRecord(BamHeader & header,
               BamIOContext<TNameStore, TNameStoreCache> & context,
               RecordReader<TStream, SinglePass<TSpec> > & reader,
               Sam const & tag)
{
    BamHeaderRecord record;
    while (nextIs(reader, SamHeader()))
    {
        clear(record);
        int res = readRecord(record, context, reader, tag);
        if (res != 0)
            return res;
        appendValue(header.records, record);

        // Get sequence information from @SQ header.
        if (record.type == BAM_HEADER_REFERENCE)
        {
            CharString sn = "unknown";
            unsigned ln = 0;
            for (unsigned i = 0; i < length(record.tags); ++i)
            {
                if (record.tags[i].i1 == "SN")
                {
                    sn = record.tags[i].i2;
                }
                else if (record.tags[i].i1 == "LN")
                {
                    if (!lexicalCast2<unsigned>(ln, record.tags[i].i2))
                        ln = 0;
                }
            }
            typedef typename BamHeader::TSequenceInfo TSequenceInfo;
            appendValue(header.sequenceInfos, TSequenceInfo(sn, ln));
            // Add name to name store cache if necessary.
            unsigned unused = 0;
            (void)unused;
            if (!getIdByName(nameStore(context), sn, unused, nameStoreCache(context)))
                appendName(nameStore(context), sn, nameStoreCache(context));
        }
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Function readRecord()                                     BamAlignmentRecord
// ----------------------------------------------------------------------------

/**
.Function.readRecord
..signature:readRecord(alignmentRecord, context, recordReader, tag)
*/

template <typename TStream, typename TSpec, typename TNameStore, typename TNameStoreCache>
int readRecord(BamAlignmentRecord & record,
               BamIOContext<TNameStore, TNameStoreCache> & context,
               RecordReader<TStream, SinglePass<TSpec> > & reader,
               Sam const & /*tag*/)
{
    clear(record);
    CharString &buffer = context.buffer;

#define SEQAN_SKIP_TAB                              \
    do                                              \
    {                                               \
        res = skipChar(reader, '\t');               \
        if (res != 0)                               \
            return res;                             \
    }                                               \
    while (false) 

    int res = 0;

    // QNAME
    res = readUntilTabOrLineBreak(record.qName, reader);
    if (res != 0)
        return res;
    SEQAN_SKIP_TAB;
    
    // FLAG
    // TODO(holtgrew): Interpret hex and char as c-samtools -X does?
    clear(buffer);
    res = readDigits(buffer, reader);
    if (res != 0)
        return res;
    record.flag = lexicalCast<__uint16>(buffer);
    SEQAN_SKIP_TAB;

    // RNAME
    clear(buffer);
    res = readUntilTabOrLineBreak(buffer, reader);
    if (res != 0)
        return res;
    if (buffer == "*")
    {
        record.rID = BamAlignmentRecord::INVALID_REFID;
    }
    else if (buffer == "0")
    {
        record.rID = BamAlignmentRecord::INVALID_REFID;
    }
    else if (!getIdByName(nameStore(context), buffer, record.rID))
    {
        record.rID = length(nameStore(context));
        appendName(nameStore(context), buffer, nameStoreCache(context));
    }
    SEQAN_SKIP_TAB;

    // POS
    clear(buffer);
    res = readUntilChar(buffer, reader, '\t');
    if (res != 0)
        return res;
    if (buffer == "*")
        record.beginPos = BamAlignmentRecord::INVALID_POS;
    else if (buffer == "0")
        record.beginPos = BamAlignmentRecord::INVALID_POS;
    else
        record.beginPos = lexicalCast<__uint32>(buffer) - 1;
    SEQAN_SKIP_TAB;

    // MAPQ
    clear(buffer);
    if (value(reader) == '*')
    {
        record.mapQ = 255;
        goNext(reader);
    }
    else
    {
        res = readDigits(buffer, reader);
        if (res != 0)
            return res;
        record.mapQ = lexicalCast<__uint16>(buffer);
    }
    SEQAN_SKIP_TAB;

    // CIGAR
    CigarElement<> element;
    if (atEnd(reader))
        return EOF_BEFORE_SUCCESS;
    if (value(reader) == '*')
    {
        goNext(reader);
    }
    else
    {
        do
        {
            clear(buffer);
            res = readDigits(buffer, reader);
            if (res != 0)
                return res;
            element.count = lexicalCast<__uint32>(buffer);
            element.operation = value(reader);
            if (goNext(reader))
                return EOF_BEFORE_SUCCESS;
            appendValue(record.cigar, element);
        } while (value(reader) != '\t');
    }
    SEQAN_SKIP_TAB;

    // RNEXT
    clear(buffer);
    res = readUntilChar(buffer, reader, '\t');
    if (res != 0)
        return res;
    if (buffer == "*")
    {
        record.rNextId = BamAlignmentRecord::INVALID_REFID;
    }
    else if (buffer == "=")
    {
        record.rNextId = record.rID;
    }
    else if (!getIdByName(nameStore(context), buffer, record.rNextId))
    {
        record.rNextId = length(nameStore(context));
        appendName(nameStore(context), buffer, nameStoreCache(context));
    }
    SEQAN_SKIP_TAB;

    // PNEXT
    if (atEnd(reader))
        return EOF_BEFORE_SUCCESS;
    if (value(reader) == '*')
    {
        record.pNext = BamAlignmentRecord::INVALID_POS;
        goNext(reader);
    }
    else
    {
        clear(buffer);
        res = readDigits(buffer, reader);
        if (res != 0)
            return res;
        if (buffer == "0")
            record.pNext = BamAlignmentRecord::INVALID_POS;
        else
            record.pNext = lexicalCast<__uint32>(buffer) - 1;
    }
    SEQAN_SKIP_TAB;

    // TLEN
    if (atEnd(reader))
        return EOF_BEFORE_SUCCESS;
    if (value(reader) == '*')
    {
        record.tLen = MaxValue<__int32>::VALUE;
        goNext(reader);
    }
    else
    {
        clear(buffer);
        if (value(reader) == '-')
        {
            appendValue(buffer, value(reader));
            if (goNext(reader))
                return SAM_INVALID_RECORD;
        }
        res = readDigits(buffer, reader);
        if (res != 0)
            return res;
        record.tLen = lexicalCast<__int32>(buffer);
    }
    SEQAN_SKIP_TAB;

    // SEQ
    res = readUntilTabOrLineBreak(record.seq, reader);
    if (res != 0)
        return res;
    // Handle case of missing sequence:  Clear seq string as documented.
    if (record.seq == "*")
        clear(record.seq);
    SEQAN_SKIP_TAB;

    // QUAL
    res = readUntilTabOrLineBreak(record.qual, reader);
    if (res == EOF_BEFORE_SUCCESS)  // The record ends on EOF.
        return 0;
    if (res != 0)
        return res;
    // Handle case of missing quality:  Clear qual string as documented.
    if (record.qual == "*")
        clear(record.qual);

    // The following list of tags is optional.  A line break or EOF could also follow.
    if (atEnd(reader))
        return 0;
    if (value(reader) != '\t')
    {
        res = skipLine(reader);
        return res;
    }
    SEQAN_SKIP_TAB;

    // TAGS
    clear(buffer);
    res = readLine(buffer, reader);
    if (res != 0 && res != EOF_BEFORE_SUCCESS)
        return res;
    assignTagsSamToBam(record.tags, buffer);

    return 0;

#undef SEQAN_SKIP_TAB
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_READ_SAM_H_
