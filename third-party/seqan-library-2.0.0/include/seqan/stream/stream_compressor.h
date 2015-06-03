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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_STREAM_STREAM_COMPRESSOR_H_
#define SEQAN_STREAM_STREAM_COMPRESSOR_H_


#if SEQAN_HAS_ZLIB
// Zlib headers
#include <zlib.h>
#include "zipstream/zutil.h"
#endif

#include <algorithm>    // copy

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TOutPager, typename TSpec>
struct Pager;

// ============================================================================
// Tags, Enums
// ============================================================================

// ============================================================================
// Classes
// ============================================================================

template <typename TAlgTag>
struct Compress;

template <typename TAlgTag>
struct CompressionContext {};

template <typename TAlgTag>
struct DefaultPageSize;

#if SEQAN_HAS_ZLIB

template <>
struct CompressionContext<GZFile>
{
    z_stream strm;

    CompressionContext()
    {
        memset(&strm, 0, sizeof(z_stream));
    }
};

template <>
struct CompressionContext<BgzfFile>:
    CompressionContext<GZFile>
{
    enum { BLOCK_HEADER_LENGTH = 18 };
    unsigned char headerPos;
};

template <typename T>
struct MagicHeader<BgzfFile, T>
{
    static char const VALUE[18];
};

template <typename T>
char const MagicHeader<BgzfFile, T>::VALUE[18] =
{
    MagicHeader<GZFile>::VALUE[0], MagicHeader<GZFile>::VALUE[1], MagicHeader<GZFile>::VALUE[2],
    4, 0, 0, 0, 0, 0, '\xff', 6, 0, 'B', 'C', 2, 0, 0, 0
};

template <>
struct DefaultPageSize<BgzfFile>
{
    static const unsigned MAX_BLOCK_SIZE = 64 * 1024;
    static const unsigned BLOCK_FOOTER_LENGTH = 8;
    // 5 bytes block overhead (see 3.2.4. at http://www.gzip.org/zlib/rfc-deflate.html)
    static const unsigned ZLIB_BLOCK_OVERHEAD = 5;

    // Reduce the maximal input size, such that the compressed data
    // always fits in one block even for level Z_NO_COMPRESSION.
    enum { BLOCK_HEADER_LENGTH = CompressionContext<BgzfFile>::BLOCK_HEADER_LENGTH };
    static const unsigned VALUE = MAX_BLOCK_SIZE - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH - ZLIB_BLOCK_OVERHEAD;
};


/*
template <typename TOutPager, typename TAlgTag>
Pager<TOutPager, Compress<TAlgTag> >
{
    TOutPager outPager;             // outbound pager
    PageTable<FixedSize> table;     // our page table

    Pager():
        table(DefaultPageSize<TAlgTag>::VALUE)
    {}

    Page & getPage (__int64 position)
    {
        Page *page;
        {
            ScopedReadLock(table.lock);

            page = table[position];
            if (posInPage(position, page))                      // does the page exist yet?
                return page;
        }
        {
            ScopedWriteLock(table.lock);

            page = table[position];
            if (posInPage(position, page))                      // does the page exist yet?
                return page;

            page = new Page(table.rangeForPos(position));       // create new page
            reserve(page.data, table.pageSize);                 // allocate required memory
            table.insertPage(page);                             // insert page
            prevPage = prevPage(position);
        }
        return page;
    }

    void putPage (Page &page)
    {
        __int64 outPosition = 0;                                // compute start position in outbound pager
        if (page.range.begin != 0)
        {
            PageRange range = getPageRange(beginPosition(page.range) - 1);
            outPosition = endPosition(range);                   // wait for end position of the previous page
        }

        TCompressionContext ctx;
        initCompressionContext(ctx);

        Size<Page>::Type leftToCompress = length(page);
        while (leftToCompress != 0)
        {
            Page &outPage = outPager.getPage(outPosition);

            auto handle = std::async(std::launch::async,
                                  parallel_sum<RAIter>, mid, end);
            compress
        }
    }
};
*/
// ============================================================================
// Functions
// ============================================================================

inline void
compressInit(CompressionContext<GZFile> & ctx)
{
    const int GZIP_WINDOW_BITS = -15;   // no zlib header
    const int Z_DEFAULT_MEM_LEVEL = 8;

    ctx.strm.zalloc = NULL;
    ctx.strm.zfree = NULL;

    // (weese:) We use Z_BEST_SPEED instead of Z_DEFAULT_COMPRESSION as it turned out
    //          to be 2x faster and produces only 7% bigger output
//    int status = deflateInit2(&ctx.strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED,
//                              GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY);
    int status = deflateInit2(&ctx.strm, Z_BEST_SPEED, Z_DEFLATED,
                              GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY);
    if (status != Z_OK)
        throw IOError("GZFile deflateInit2() failed.");
}

inline void
compressInit(CompressionContext<BgzfFile> & ctx)
{
    compressInit(static_cast<CompressionContext<GZFile> &>(ctx));
    ctx.headerPos = 0;
}

template <typename TTarget, typename TSourceIterator>
inline typename Size<TTarget>::Type
compress(TTarget & target, TSourceIterator & source, CompressionContext<BgzfFile> & ctx)
{
    typedef typename Chunk<TTarget>::Type           TTargetChunk;
    typedef typename Chunk<TSourceIterator>::Type   TSourceChunk;
    typedef typename Value<TSourceChunk>::Type      TSourceValue;

    TTargetChunk tChunk;
    TSourceChunk sChunk;

    if (ctx.headerPos < sizeof(MagicHeader<BgzfFile>::VALUE))
    {
        size_t headerLeft = sizeof(MagicHeader<BgzfFile>::VALUE) - ctx.headerPos;
        reserveChunk(target, headerLeft, Output());

        tChunk = getChunk(target, Output());
        size_t size = std::min(headerLeft, length(tChunk));
        SEQAN_ASSERT_GT(size, 0u);

        std::copy(tChunk.begin, sChunk.begin, size);

        advanceChunk(target, size);
        ctx.headerPos += size;
        return size;
    }
    else
    {
        sChunk = getChunk(source, Input());
        tChunk = getChunk(target, Output());

        ctx.strm.next_in = static_cast<Bytef *>(sChunk.begin);
        ctx.strm.next_out = static_cast<Bytef *>(tChunk.begin);
        ctx.strm.avail_in = length(sChunk) * sizeof(TSourceValue);
        ctx.strm.avail_out = length(tChunk);

        SEQAN_ASSERT_GT(ctx.strm.avail_out, 0u);

        int status = deflate(&ctx.strm, Z_NO_FLUSH);
        if (status != Z_OK)
            throw IOError("BgzfFile deflateInit2() failed.");

        source += length(sChunk) - ctx.strm.avail_in;
        size_t size = length(tChunk) - ctx.strm.avail_out;
        advanceChunk(target, size);
        return size;
    }


//    status = deflate(&zs, Z_FINISH);
//    bool rawDataTooBig = (status != Z_STREAM_END);
//
//    status = deflateEnd(&zs);
//    if (status != Z_OK)
//        throw IOError("BgzfFile deflateEnd() failed.");
//
//    if (!rawDataTooBig)
//    {
//        resize(page.raw, zs.total_out + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH);
//        break;
//    }
}

// ----------------------------------------------------------------------------
// Helper Function _bgzfUnpackXX()
// ----------------------------------------------------------------------------

inline unsigned short
_bgzfUnpack16(char const * buffer)
{
    return *reinterpret_cast<unsigned short const *>(buffer);
}

inline unsigned
_bgzfUnpack32(char const * buffer)
{
    return *reinterpret_cast<unsigned const *>(buffer);
}

// ----------------------------------------------------------------------------
// Helper Function _bgzfPackXX()
// ----------------------------------------------------------------------------

inline void
_bgzfPack16(char * buffer, unsigned short value)
{
    *reinterpret_cast<unsigned short *>(buffer) = value;
}

inline void
_bgzfPack32(char * buffer, unsigned value)
{
    *reinterpret_cast<unsigned *>(buffer) = value;
}


template <typename TDestValue, typename TDestCapacity, typename TSourceValue, typename TSourceLength>
inline TDestCapacity
_compressBlock(TDestValue *dstBegin,   TDestCapacity dstCapacity,
               TSourceValue *srcBegin, TSourceLength srcLength, CompressionContext<BgzfFile> & ctx)
{
    const size_t BLOCK_HEADER_LENGTH = DefaultPageSize<BgzfFile>::BLOCK_HEADER_LENGTH;
    const size_t BLOCK_FOOTER_LENGTH = DefaultPageSize<BgzfFile>::BLOCK_FOOTER_LENGTH;

    SEQAN_ASSERT_GT(dstCapacity, BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH);
    SEQAN_ASSERT_EQ(sizeof(TDestValue), 1u);
    SEQAN_ASSERT_EQ(sizeof(unsigned), 4u);

    // 1. COPY HEADER

    std::copy(&MagicHeader<BgzfFile>::VALUE[0], &MagicHeader<BgzfFile>::VALUE[BLOCK_HEADER_LENGTH], dstBegin);


    // 2. COMPRESS

    compressInit(ctx);
    ctx.strm.next_in = (Bytef *)(srcBegin);
    ctx.strm.next_out = (Bytef *)(dstBegin + BLOCK_HEADER_LENGTH);
    ctx.strm.avail_in = srcLength * sizeof(TSourceValue);
    ctx.strm.avail_out = dstCapacity - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;

    int status = deflate(&ctx.strm, Z_FINISH);
    if (status != Z_STREAM_END)
    {
        deflateEnd(&ctx.strm);
        throw IOError("Deflation failed. Compressed BGZF data is too big.");
    }

    status = deflateEnd(&ctx.strm);
    if (status != Z_OK)
        throw IOError("BGZF deflateEnd() failed.");


    // 3. APPEND FOOTER

    // Set compressed length into buffer, compute CRC and write CRC into buffer.

    size_t len = dstCapacity - ctx.strm.avail_out;
    _bgzfPack16(dstBegin + 16, len - 1);

    dstBegin += len - BLOCK_FOOTER_LENGTH;
    _bgzfPack32(dstBegin, crc32(crc32(0u, NULL, 0u), (Bytef *)(srcBegin), srcLength * sizeof(TSourceValue)));
    _bgzfPack32(dstBegin + 4, srcLength * sizeof(TSourceValue));

    return dstCapacity - ctx.strm.avail_out;
}

inline void
decompressInit(CompressionContext<GZFile> & ctx)
{
    const int GZIP_WINDOW_BITS = -15;   // no zlib header

    ctx.strm.zalloc = NULL;
    ctx.strm.zfree = NULL;
    int status = inflateInit2(&ctx.strm, GZIP_WINDOW_BITS);
    if (status != Z_OK)
        throw IOError("GZip inflateInit2() failed.");
}

inline void
decompressInit(CompressionContext<BgzfFile> & ctx)
{
    decompressInit(static_cast<CompressionContext<GZFile> &>(ctx));
    ctx.headerPos = 0;
}

inline bool
_bgzfCheckHeader(char const * header)
{
    const char FLG_FEXTRA = 4;
    const char BGZF_ID1 = 'B';
    const char BGZF_ID2 = 'C';
    const char BGZF_LEN = 2;
    const char BGZF_XLEN = 6;  // BGZF_LEN+4

    return (header[0] == (char)MagicHeader<GZFile>::VALUE[0] &&
            header[1] == (char)MagicHeader<GZFile>::VALUE[1] &&
            header[2] == (char)MagicHeader<GZFile>::VALUE[2] &&
            (header[3] & FLG_FEXTRA) != 0 &&
            _bgzfUnpack16(header + 10) == BGZF_XLEN &&
            header[12] == BGZF_ID1 &&
            header[13] == BGZF_ID2 &&
            _bgzfUnpack16(header + 14) == BGZF_LEN);
}

// read first bytes of a file/stream and compare with file format's magic header
template <typename TStream>
inline bool
guessFormatFromStream(TStream &istream, BgzfFile)
{
    char putbackBuf[18];
    bool match = false;

    SEQAN_ASSERT(istream.good());

    // try to read and check header
    size_t numRead = istream.readsome(&putbackBuf[0], sizeof(putbackBuf));
    if (numRead == sizeof(putbackBuf) && _bgzfCheckHeader(putbackBuf))
        match = true;

    // unget all read characters
    for (; numRead > 0; --numRead)
        istream.unget();

    SEQAN_ASSERT(istream.good());

    return match;
}

// ----------------------------------------------------------------------------
// Function _preprocessFilePage()
// ----------------------------------------------------------------------------

template <typename TDestValue, typename TDestCapacity, typename TSourceValue, typename TSourceLength>
inline TDestCapacity
_decompressBlock(TDestValue *dstBegin,   TDestCapacity dstCapacity,
                 TSourceValue *srcBegin, TSourceLength srcLength, CompressionContext<BgzfFile> & ctx)
{
    const size_t BLOCK_HEADER_LENGTH = DefaultPageSize<BgzfFile>::BLOCK_HEADER_LENGTH;
    const size_t BLOCK_FOOTER_LENGTH = DefaultPageSize<BgzfFile>::BLOCK_FOOTER_LENGTH;

    SEQAN_ASSERT_EQ(sizeof(TSourceValue), 1u);
    SEQAN_ASSERT_EQ(sizeof(unsigned), 4u);

    // 1. CHECK HEADER

    if (srcLength <= BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH)
        throw IOError("BGZF block too short.");

    if (!_bgzfCheckHeader(srcBegin))
        throw IOError("Invalid BGZF block header.");

    size_t compressedLen = _bgzfUnpack16(srcBegin + 16) + 1u;
    if (compressedLen != srcLength)
        throw IOError("BGZF compressed size mismatch.");


    // 2. DECOMPRESS

    decompressInit(ctx);
    ctx.strm.next_in = (Bytef *)(srcBegin + BLOCK_HEADER_LENGTH);
    ctx.strm.next_out = (Bytef *)(dstBegin);
    ctx.strm.avail_in = srcLength - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;
    ctx.strm.avail_out = dstCapacity * sizeof(TDestValue);

    int status = inflate(&ctx.strm, Z_FINISH);
    if (status != Z_STREAM_END)
    {
        inflateEnd(&ctx.strm);
        throw IOError("Inflation failed. Decompressed BGZF data is too big.");
    }

    status = inflateEnd(&ctx.strm);
    if (status != Z_OK)
        throw IOError("BGZF inflateEnd() failed.");


    // 3. CHECK FOOTER

    // Check compressed length in buffer, compute CRC and compare with CRC in buffer.

    unsigned crc = crc32(crc32(0u, NULL, 0u), (Bytef *)(dstBegin), dstCapacity - ctx.strm.avail_out);

    srcBegin += compressedLen - BLOCK_FOOTER_LENGTH;
    if (_bgzfUnpack32(srcBegin) != crc)
        throw IOError("BGZF wrong checksum.");

    if (_bgzfUnpack32(srcBegin + 4) != dstCapacity - ctx.strm.avail_out)
        throw IOError("BGZF size mismatch.");

    return (dstCapacity - ctx.strm.avail_out) / sizeof(TDestValue);
}

#endif  // #if SEQAN_HAS_ZLIB

}  // namespace seqan

#endif  // SEQAN_STREAM_STREAM_COMPRESSOR_H_
