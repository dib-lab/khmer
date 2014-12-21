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
// Code for writing BAM.
// ==========================================================================

// TODO(holtgrew): Add buffer to context?

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_WRITE_BAM_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_WRITE_BAM_H_

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
// Function writeRecord()                                             BamHeader
// ----------------------------------------------------------------------------

template <typename TStream, typename TNameStore, typename TNameStoreCache>
int write2(TStream & stream,
           BamHeader const & header,
           BamIOContext<TNameStore, TNameStoreCache> const & context,
           Bam const & /*tag*/)
{
    int res = streamWriteBlock(stream, "BAM\1", 4);
    if (res != 4)
        return 1;  // Could not write magic.

    // Create text of header.
    CharString headerBuffer;
    for (unsigned i = 0; i < length(header.records); ++i)
    {
        res = write2(headerBuffer, header.records[i], context, Sam());
        if (res != 0)
            return 1;  // Error writing header to buffer.
    }
    // Note that we do not write out a null-character to terminate the header.  This would be valid by the SAM standard
    // but the samtools do not expect this and write out the '\0' when converting from BAM to SAM.
    // appendValue(headerBuffer, '\0');

    // Write text header.
    __int32 lText = length(headerBuffer);
    res = streamWriteBlock(stream, reinterpret_cast<char const *>(&lText), 4);
    if (res != 4)
        return 1;  // Error writing l_text.

    res = streamWriteBlock(stream, &headerBuffer[0], lText);

    // Write references.
    __int32 nRef = _max(length(header.sequenceInfos), length(nameStore(context)));
    res = streamWriteBlock(stream, reinterpret_cast<char const *>(&nRef), 4);
    if (res != 4)
        return 1;  // Error writing n_ref;

    for (unsigned i = 0; i < length(header.sequenceInfos); ++i)
    {
        __int32 lName = length(header.sequenceInfos[i].i1) + 1;
        res = streamWriteBlock(stream, reinterpret_cast<char const *>(&lName), 4);
        if (res != 4)
            return 1;  // Error writing l_name;

        res = streamWriteBlock(stream, &header.sequenceInfos[i].i1[0], lName - 1);
        if (res != lName - 1)
            return 1;  // Error writing name;

        char const n = '\0';
        res = streamWriteBlock(stream, &n, 1);
        if (res != 1)
            return 1;  // Error writing trailing '\0'.

        __int32 lRef = header.sequenceInfos[i].i2;
        res = streamWriteBlock(stream, reinterpret_cast<char const *>(&lRef), 4);
        if (res != 4)
            return 1;  // Error writing l_ref;
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Function writeRecord()                                    BamAlignmentRecord
// ----------------------------------------------------------------------------

static inline int _reg2Bin(uint32_t beg, uint32_t end)
{
    --end;
    if (beg >> 14 == end >> 14)
        return 4681 + (beg >> 14);

    if (beg >> 17 == end >> 17)
        return 585 + (beg >> 17);

    if (beg >> 20 == end >> 20)
        return 73 + (beg >> 20);

    if (beg >> 23 == end >> 23)
        return 9 + (beg >> 23);

    if (beg >> 26 == end >> 26)
        return 1 + (beg >> 26);

    return 0;
}

template <typename TStream, typename TNameStore, typename TNameStoreCache>
int write2(TStream & stream,
           BamAlignmentRecord const & record,
           BamIOContext<TNameStore, TNameStoreCache> const & /*context*/,
           Bam const & /*tag*/)
{
    CharString buffer;

    // First, write record to buffer.

    // refID
    streamWriteBlock(buffer, reinterpret_cast<char const *>(&record.rID), 4);

    // pos
    streamWriteBlock(buffer, reinterpret_cast<char const *>(&record.beginPos), 4);

    // bin_mq_nl
    SEQAN_ASSERT_LT(length(record.qName) + 1u, 255u);
    __uint8 lReadName = length(record.qName) + 1;
    unsigned l = 0;
    _getLengthInRef(record.cigar, l);
    __uint32 bin = _reg2Bin(record.beginPos, record.beginPos + l);
    __uint32 binMqNl = (bin << 16) | (record.mapQ << 8) | lReadName;
    streamWriteBlock(buffer, reinterpret_cast<char const *>(&binMqNl), 4);

    // flag_nc
    __uint16 nCigarOp = length(record.cigar);
    __uint32 flagNc = (record.flag << 16) | nCigarOp;
    streamWriteBlock(buffer, reinterpret_cast<char const *>(&flagNc), 4);

    // l_seq
    __int32 lSeq = length(record.seq);
    streamWriteBlock(buffer, reinterpret_cast<char const *>(&lSeq), 4);

    // next_refID
    streamWriteBlock(buffer, reinterpret_cast<char const *>(&record.rNextId), 4);

    // next_pos
    streamWriteBlock(buffer, reinterpret_cast<char const *>(&record.pNext), 4);

    // tlen
    __int32 zero = 0;
    if (record.tLen == BamAlignmentRecord::INVALID_LEN)
        streamWriteBlock(buffer, reinterpret_cast<char const *>(&zero), 4);
    else
        streamWriteBlock(buffer, reinterpret_cast<char const *>(&record.tLen), 4);

    // read_name
    streamWriteBlock(buffer, reinterpret_cast<char const *>(&record.qName[0]), lReadName - 1);
    streamWriteChar(buffer, '\0');

    // cigar
    static __uint8 const MAP[256] =
    {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0,
        0, 0, 0, 0, 2, 0, 0, 0, 5, 1, 0, 0, 0, 0, 3, 0,
        6, 0, 0, 4, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };
    for (unsigned i = 0; i < length(record.cigar); ++i)
    {
        __uint32 x = record.cigar[i].count;
        x <<= 4;
        x |= MAP[static_cast<int>(record.cigar[i].operation)];
        streamWriteBlock(buffer, reinterpret_cast<char const *>(&x), 4);
    }

    // seq
    static __uint8 const MAP2[256] =
    {
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 0, 15, 15,
        15, 1, 14, 2, 13, 15, 15, 4, 11, 15, 15, 12, 15, 3, 15, 15,
        15, 15, 5, 6, 8, 15, 7, 9, 15, 10, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15
    };
    __uint8 c = 0;
    for (int i = 0; i < lSeq; ++i)
    {
        c <<= 4;
        c &= 0xf0;
        c |= MAP2[static_cast<int>(record.seq[i])];
        if (i % 2 == 1)
            streamWriteChar(buffer, c);
    }
    if (lSeq % 2 == 1)
    {
        c <<= 4;
        c &= 0xf0;
        streamWriteChar(buffer, c);
    }

    // qual
    if (empty(record.qual))
    {
        for (unsigned i = 0; i < length(record.qual); ++i)
            streamWriteChar(buffer, static_cast<unsigned char>(0xff));
    }
    else
    {
        for (unsigned i = 0; i < length(record.qual); ++i)
            streamWriteChar(buffer, static_cast<char>(record.qual[i] - '!'));
    }

    // tags
    if (length(record.tags) > 0u)
        streamWriteBlock(buffer, reinterpret_cast<char const *>(&record.tags[0]), length(record.tags));

    // buffer to stream
    __uint32 blockSize = length(buffer);
    streamWriteBlock(stream, reinterpret_cast<char const *>(&blockSize), 4);
    return streamWriteBlock(stream, &buffer[0], blockSize) != blockSize;
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_WRITE_BAM_H_
