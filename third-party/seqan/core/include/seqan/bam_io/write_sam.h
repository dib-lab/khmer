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
// Code for writing SAM.
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_WRITE_SAM_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_WRITE_SAM_H_

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

// ----------------------------------------------------------------------------
// Function write2()                                            BamHeaderRecord
// ----------------------------------------------------------------------------

template <typename TStream, typename TNameStore, typename TNameStoreCache>
int write2(TStream & stream,
           BamHeaderRecord const & header,
           BamIOContext<TNameStore, TNameStoreCache> const & /*context*/,
           Sam const & /*tag*/)
{
    char const * headerTypes[] = {"@HD", "@SQ", "@RG", "@PG", "@CO"};
    streamPut(stream, headerTypes[header.type]);
    if (header.type == BAM_HEADER_COMMENT)
    {
        streamPut(stream, header.tags[0].i2);
    }
    else
    {
        for (unsigned i = 0; i < length(header.tags); ++i)
        {
            streamPut(stream, '\t');
            streamPut(stream, header.tags[i].i1);
            streamPut(stream, ':');
            streamPut(stream, header.tags[i].i2);
        }
    }

    int res = streamPut(stream, '\n');
    if (res != 0)
        return res;

    return 0;
}

// ----------------------------------------------------------------------------
// Function write2()                                                  BamHeader
// ----------------------------------------------------------------------------

template <typename TStream, typename TNameStore, typename TNameStoreCache>
int write2(TStream & stream,
           BamHeader const & header,
           BamIOContext<TNameStore, TNameStoreCache> const & context,
           Sam const & tag)
{
    std::set<CharString> writtenSeqInfos;

    for (unsigned i = 0; i < length(header.records); ++i)
    {
        BamHeaderRecord const & record = header.records[i];
        if (record.type == BAM_HEADER_REFERENCE)
        {
            for (unsigned i = 0; i < length(record.tags); ++i)
            {
                if (record.tags[i].i1 == "SN")
                {
                    writtenSeqInfos.insert(record.tags[i].i2);
                    break;
                }
            }
        }

        int res = write2(stream, record, context, tag);
        if (res != 0)
            return res;
    }

    // Write missing @SQ header records.
    for (unsigned i = 0; i < length(header.sequenceInfos); ++i)
    {
        if (writtenSeqInfos.find(header.sequenceInfos[i].i1) != writtenSeqInfos.end())
            continue;
        int res = streamPut(stream, "@SQ\tSN:");
        if (res != 0)
            return res;

        res = streamPut(stream, header.sequenceInfos[i].i1);
        if (res != 0)
            return res;

        res = streamPut(stream, "\tLN:");
        if (res != 0)
            return res;

        res = streamPut(stream, header.sequenceInfos[i].i2);
        if (res != 0)
            return res;

        res = streamPut(stream, '\n');
        if (res != 0)
            return res;
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Function write2()                                         BamAlignmentRecord
// ----------------------------------------------------------------------------

template <typename TStream, typename TNameStore, typename TNameStoreCache>
int write2(TStream & stream,
           BamAlignmentRecord const & record,
           BamIOContext<TNameStore, TNameStoreCache> const & context,
           Sam const & /*tag*/)
{
    int res = 0;

#define SEQAN_PUT_TAB                           \
    do {                                        \
        res = streamPut(stream, '\t');      \
        if (res != 0)                       \
            return res;                     \
    } \
    while (false)

    res = streamPut(stream, record.qName);
    if (res != 0)
        return res;

    SEQAN_PUT_TAB;

    res = streamPut(stream, record.flag);
    if (res != 0)
        return res;

    SEQAN_PUT_TAB;

    if (record.rID == BamAlignmentRecord::INVALID_REFID)
        res = streamPut(stream, '*');
    else
        res = streamPut(stream, nameStore(context)[record.rID]);
    if (res != 0)
        return res;

    SEQAN_PUT_TAB;

    if (record.rID == BamAlignmentRecord::INVALID_REFID)
        res = streamPut(stream, '*');
    else
        res = streamPut(stream, record.beginPos + 1);
    if (res != 0)
        return res;

    SEQAN_PUT_TAB;

    res = streamPut(stream, static_cast<__uint16>(record.mapQ));
    if (res != 0)
        return res;

    SEQAN_PUT_TAB;

    if (length(record.cigar) == 0u)
    {
        res = streamPut(stream, '*');
        if (res != 0)
            return res;
    }
    else
    {
        for (unsigned i = 0; i < length(record.cigar); ++i)
        {
            res = streamPut(stream, record.cigar[i].count);
            if (res != 0)
                return res;

            res = streamPut(stream, record.cigar[i].operation);
            if (res != 0)
                return res;
        }
    }

    SEQAN_PUT_TAB;

    if (record.rNextId == BamAlignmentRecord::INVALID_REFID)
        res = streamPut(stream, '*');
    else if (record.rID == record.rNextId)
        res = streamPut(stream, '=');
    else
        res = streamPut(stream, nameStore(context)[record.rNextId]);
    if (res != 0)
        return res;

    SEQAN_PUT_TAB;

    if (record.pNext == BamAlignmentRecord::INVALID_POS)
        res = streamPut(stream, '0');
    else
        res = streamPut(stream, record.pNext + 1);
    if (res != 0)
        return res;

    SEQAN_PUT_TAB;

    if (record.tLen == BamAlignmentRecord::INVALID_LEN)
        res = streamPut(stream, '0');
    else
        res = streamPut(stream, record.tLen);
    if (res != 0)
        return res;

    SEQAN_PUT_TAB;

    if (empty(record.seq))
        res = streamPut(stream, '*');  // Case of empty seq string / "*".
    else
        res = streamPut(stream, record.seq);
    if (res != 0)
        return res;

    SEQAN_PUT_TAB;


    if (empty(record.qual))  // Case of empty quality string / "*".
        res = streamPut(stream, '*');
    else
        res = streamPut(stream, record.qual);
    if (res != 0)
        return res;

    if (length(record.tags) > 0u)
    {
        SEQAN_PUT_TAB;
        CharString buffer;
        assignTagsBamToSam(buffer, record.tags);
        streamPut(stream, buffer);
    }

    return streamPut(stream, '\n');

#undef SEQAN_PUT_TAB
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_WRITE_SAM_H_
