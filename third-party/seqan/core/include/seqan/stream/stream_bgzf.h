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
// BGZF Stream Access.
//
// This code is based on the original code by the Broad institute at MIT.
// The license for the original code is reproduced below.
// ==========================================================================

/* The MIT License

   Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#ifndef EXTRAS_INCLUDE_SEQAN_STREAM_STREAM_BGZF_H_
#define EXTRAS_INCLUDE_SEQAN_STREAM_STREAM_BGZF_H_

#include <map>

#include <zlib.h>

// TODO(holtgrew): Check file for truncation via function.

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

struct Bgzf_;
typedef Tag<Bgzf_> Bgzf;
template <> class Stream<Bgzf>;
inline void close(Stream<Bgzf> & stream);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// One entry in the BGZF cache.

struct BgzfCacheEntry_
{
    int size;
    String<char> block;
    __int64 endOffset;
};

/**
.Spec.BGZF Stream
..cat:Input/Output
..signature:Stream<Bgzf>
..general:Class.Stream
..summary:Reading and writing of BGZF data.
..remarks:Not copy constructable, not assignable.
..remarks:
BGZF is the Block GZip Format which is used as the underlying format for BAM and TABIX.
Data is written out compressed with gzip but the uncompressed data is split into blocks with a maximum block size.
It is therefore possible to jump to beginnings of blocks in the resulting files, decompress the block and then jump into the block itself.
..include:seqan/stream.h
..example.code:
Stream<Bgzf> stream;
if (!open(stream, "myfile.bam", "r"))
    return 1;  // error
// ... work
 */

struct Bgzf_;
typedef Tag<Bgzf_> Bgzf;

template <>
class Stream<Bgzf>
{
public:
    typedef File<Sync<> > TFile;  // Type to use for reading/writing files.
    typedef Position<TFile>::Type TPosition;

    // The current error state, 0 for no error.
    int _error;

    // Flag that indicates whether we are at the end of the file.
    bool _atEof;

    // Mode of currently opened file.
    int _openMode;  // TODO(holtgrew): Should go into File class.

    // Current zlib compression level.
    int _compressLevel;

    // The file handle for the BGZF stream.
    TFile _file;

    // Buffer for uncompressed data.
    String<char> _uncompressedBlock;

    // Buffer for compressed data.
    String<char> _compressedBlock;

    // The address of the current block.
    TPosition _blockPosition;

    // Length of the current block.
    __int32 _blockLength;

    // Offset in the current block.
    __int32 _blockOffset;

    // Cache of decompressed blocks.
    std::map<__int64, BgzfCacheEntry_ *> _cache;

    // Number of bytes in cached blocks.
    int _cacheSize;

    // Maximum cache size, as number of bytes in cached blocks.
    int _maxCacheSize;

    // Whether or not the file is owned (i.e. opened with open()) or just attached to an already open file via POSIX
    // file handle.
    bool _fileOwned;

    // Size of the file in bytes as it is on the disk.
    __int64 _fileSize;

    Stream() : _error(0), _atEof(false), _openMode(0), _compressLevel(Z_DEFAULT_COMPRESSION), _blockPosition(0),
               _blockLength(0), _blockOffset(0), _cacheSize(0), _maxCacheSize(0), _fileOwned(false), _fileSize(0)
    {}

    ~Stream()
    {
        close(*this);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

template <>
struct Difference<Stream<Bgzf> >
{
    typedef __int64 Type;
};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <>
struct Position<Stream<Bgzf> >
{
    typedef __int64 Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <>
struct Size<Stream<Bgzf> >
{
    typedef __int64 Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <>
struct Value<Stream<Bgzf> >
{
    typedef char Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, IsInput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<Bgzf>, IsInput>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, IsOutput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<Bgzf>, IsOutput>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, HasPeek>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<Bgzf>, HasPeek>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, HasFilename>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<Bgzf>, HasFilename>
{
    typedef False Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, Seek<TSpec> >
// ----------------------------------------------------------------------------

template <typename TSpec>
struct HasStreamFeature<Stream<Bgzf>, Seek<TSpec> >
{
    typedef False Type;
};

template <>
struct HasStreamFeature<Stream<Bgzf>, Seek<OriginBegin> >
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, Tell>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<Bgzf>, Tell>
{
    typedef True Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function size()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Do we want a trait for this? Only possible for random access streams and those that have an underlying file.
///.Function.size.param.file.type:Spec.BGZF Stream
///.Function.size.class:Spec.BGZF Stream

// Just forward to underlying Class.File's size().

/*
inline typename Size<Stream<Bgzf> >::Type
size(Stream<Bgzf> const & stream)
{
    return size(stream._file);
}
*/

// ----------------------------------------------------------------------------
// Helper Function _bgzfUnpackInt16()
// ----------------------------------------------------------------------------

inline int
_bgzfUnpackInt16(unsigned char const * buffer)
{
    return (buffer[0] | (buffer[1] << 8));
}

// ----------------------------------------------------------------------------
// Helper Function _bgzfPackInt16()
// ----------------------------------------------------------------------------

inline
void
_bgzfPackInt16(unsigned char * buffer, __uint16 value)
{
    buffer[0] = value;
    buffer[1] = value >> 8;
}

// ----------------------------------------------------------------------------
// Helper Function _bgzfPackInt32()
// ----------------------------------------------------------------------------

inline void
_bgzfPackInt32(unsigned char * buffer, __int32 value)
{
    buffer[0] = value;
    buffer[1] = value >> 8;
    buffer[2] = value >> 16;
    buffer[3] = value >> 24;
}

// ----------------------------------------------------------------------------
// Helper Function _bgzfCheckHeader()
// ----------------------------------------------------------------------------

inline int
_bgzfCheckHeader(char const * header)
{
    const char GZIP_ID1 = 31;
    const __uint8 GZIP_ID2 = 139;
    const char FLG_FEXTRA = 4;
    const char BGZF_ID1 = 66;  // 'B'
    const char BGZF_ID2 = 67;  // 'C'
    const char BGZF_LEN = 2;
    const char BGZF_XLEN = 6;  // BGZF_LEN+4

    return (header[0] == GZIP_ID1 &&
            (__uint8)header[1] == (__uint8) GZIP_ID2 &&
            header[2] == Z_DEFLATED &&
            (header[3] & FLG_FEXTRA) != 0 &&
            _bgzfUnpackInt16((unsigned char*)&header[10]) == BGZF_XLEN &&
            header[12] == BGZF_ID1 &&
            header[13] == BGZF_ID2 &&
            _bgzfUnpackInt16((unsigned char*)&header[14]) == BGZF_LEN);
}

// ----------------------------------------------------------------------------
// Helper Function _bgzfLoadBlockFromCache()
// ----------------------------------------------------------------------------

// Returns number of bytes loaded from cache, 0 if address could not be found.

inline int
_bgzfLoadBlockFromCache(Stream<Bgzf> & stream, __int64 blockAddress)
{
    // If there is no block in the cache with this address then return false.
    std::map<__int64, BgzfCacheEntry_ *>::iterator it = stream._cache.find(blockAddress);
    if (it == stream._cache.end())
        return 0;

    // Update fields of stream.
    if (stream._blockLength != 0)
        stream._blockOffset = 0;
    stream._blockPosition = blockAddress;
    stream._blockLength = it->second->size;

    // Copy data from cache into uncompressed block buffer.
    if (!empty(it->second->block))
        memcpy(&stream._uncompressedBlock[0], &it->second->block[0], length(it->second->block));

    // Seek to end of cached block in the underlying file.
    seek(stream._file, it->second->endOffset, SEEK_SET);

    return it->second->size;
}

// ----------------------------------------------------------------------------
// Helper Function _bgzfCacheBlock()
// ----------------------------------------------------------------------------

inline bool
_bgzfCacheBlock(Stream<Bgzf> & stream, size_t size)
{
    if (stream._blockLength > stream._maxCacheSize)
        return false;  // Cannot cache this block.
    if (stream._cache.find(stream._blockPosition) != stream._cache.end())
        return true;  // Block is already cached.

    // TODO(holtgrew): Implement LRU scheme?
    
    // Throw out blocks from cache until the new cache fits.  We simply remove the first block which should not make
    // performance degrade much in practice.
    while (stream._cacheSize + stream._blockLength > stream._maxCacheSize)
    {
        SEQAN_ASSERT_NOT(stream._cache.empty());
        std::map<__int64, BgzfCacheEntry_ *>::iterator it = stream._cache.begin();
        stream._cacheSize -= length(it->second->block);
        delete it->second;
        stream._cache.erase(it);
    }

    // Create new cache entry, copy out block data and put it into the cache.
    BgzfCacheEntry_ * entry = new BgzfCacheEntry_();
    entry->size = stream._blockLength;
    entry->block = prefix(stream._uncompressedBlock, stream._blockLength);
    entry->endOffset = stream._blockPosition + size;
    stream._cache[stream._blockPosition] = entry;

    stream._cacheSize += length(entry->block);
    SEQAN_ASSERT_LEQ(stream._cacheSize, stream._maxCacheSize);

    return 0;
}

// ----------------------------------------------------------------------------
// Helper Function _bgzfCacheBlock()
// ----------------------------------------------------------------------------

inline void
_bgzfClearCache(Stream<Bgzf> & stream)
{
    typedef std::map<__int64, BgzfCacheEntry_ *>::iterator TIterator;
    for (TIterator it = stream._cache.begin(); it != stream._cache.end(); ++it)
        delete it->second;
    stream._cacheSize = 0;
    stream._cache.clear();
}

// ----------------------------------------------------------------------------
// Helper Function _bgzfInflateBlock()
// ----------------------------------------------------------------------------

// Inflate from compression to decompression buffer.

inline int
_bgzfInflateBlock(Stream<Bgzf> & stream, size_t blockLength)
{
    int const GZIP_WINDOW_BITS = -15;  // no zlib header

    z_stream zs;
	int status;
    zs.zalloc = NULL;
    zs.zfree = NULL;
    zs.next_in = static_cast<Bytef *>(static_cast<void *>(&stream._compressedBlock[0])) + 18;
    zs.avail_in = blockLength - 16;
    zs.next_out = static_cast<Bytef *>(static_cast<void *>(&stream._uncompressedBlock[0]));
    zs.avail_out = length(stream._uncompressedBlock);

    status = inflateInit2(&zs, GZIP_WINDOW_BITS);
    if (status != Z_OK)
        return -1;  // inflateInit2() has failed.

    status = inflate(&zs, Z_FINISH);
    if (status != Z_STREAM_END)
    {
        inflateEnd(&zs);
        return -1;  // inflate() has failed.
    }

    status = inflateEnd(&zs);
    if (status != Z_OK)
        return -1;  // inflateEnd() failed.

    return zs.total_out;
}

// ----------------------------------------------------------------------------
// Helper Function _bgzfReadBlock()
// ----------------------------------------------------------------------------

// Returns 0 on success, -1 on error, -2 on eof.

inline int
_bgzfReadBlock(Stream<Bgzf> & stream)
{
    int const BLOCK_HEADER_LENGTH = 18;
    // Make sure there is enough space in the buffer for compressed data.
    unsigned const MAX_BLOCK_SIZE = 64 * 1024;
    resize(stream._compressedBlock, MAX_BLOCK_SIZE);
    resize(stream._uncompressedBlock, MAX_BLOCK_SIZE);

    char header[BLOCK_HEADER_LENGTH];

    // Get address from block and try to get cached block from this address.
    __int64 blockAddress = tell(stream._file);
    if (_bgzfLoadBlockFromCache(stream, blockAddress))
        return 0;

    // Try to read the heder.
    __int64 posBefore = tell(stream._file);
    // TODO(holtgrew): Complicated reading because File<> interface is not so good.
    bool success = read(stream._file, &header[0], sizeof(header));
    int count = tell(stream._file) - posBefore;
    if (!success && count == 0)
        return -2;  // EOF.
    if (!success)  // Unnecessary, could go away when file interface changes.
        return -1; // Could not read header.

    // If no data could be read for the header then we are at the end of the file, this is no error.
    // TODO(holtgrew): Correct with EOF?
    if (count == 0)
    {
        stream._blockLength = 0;
        return 0;
    }

    int size = count;

    // Check that the header is valid.
    if (count != sizeof(header))
        return -1;  // Could not read the full header.
    if (!_bgzfCheckHeader(header))
        return -1;  // Header was invalid.

    // Copy header into buffer for compressed data.
    int blockLength = _bgzfUnpackInt16((unsigned char *)&header[16]) + 1;
    char * compressedBlock = (char *)&stream._compressedBlock[0];
    memcpy(compressedBlock, header, BLOCK_HEADER_LENGTH);
    int remaining = blockLength - BLOCK_HEADER_LENGTH;

    // Read remainder of block into buffer for compressed data.
    // TODO(holtgrew): Complicated reading because File<> interface is not so good.
    posBefore = tell(stream._file);
    success = read(stream._file, &compressedBlock[BLOCK_HEADER_LENGTH], remaining);
    if (!success)
        return -1;
    count = tell(stream._file) - posBefore;
    if (count != remaining)
        return -1;  // Read failed.

    size += count;

    // Decompress between compression and decompression buffer.
    count = _bgzfInflateBlock(stream, blockLength);
    if (count < 0)
        return -1;  // Decompression failed.

    if (stream._blockLength != 0)
        stream._blockOffset = 0;  // Do not reset offset if this read follows a seek.

    // Update block address and length in stream object and add block to cache.
    stream._blockPosition = blockAddress;
    stream._blockLength = count;
    _bgzfCacheBlock(stream, size);

    return 0;
}

// ----------------------------------------------------------------------------
// Helper Function _bgzfDeflateBlock()
// ----------------------------------------------------------------------------

// Deflate from uncompressed block to compressed block.  Also add extra field that stores the compressed block length.

inline int
_bgzfDeflateBlock(Stream<Bgzf> & stream, int blockLength)
{
    const int BLOCK_HEADER_LENGTH = 18;
    const int BLOCK_FOOTER_LENGTH = 8;

    const char GZIP_ID1 = 31;
    const __uint8 GZIP_ID2 = 139;
    const char CM_DEFLATE = 8;
    const char FLG_FEXTRA = 4;
    const __uint8 OS_UNKNOWN = 255;
    const char BGZF_ID1 = 66;  // 'B'
    const char BGZF_ID2 = 67;  // 'C'
    const char BGZF_LEN = 2;
    const char BGZF_XLEN = 6;  // BGZF_LEN+4

    const int GZIP_WINDOW_BITS = -15; // no zlib header
    const int Z_DEFAULT_MEM_LEVEL = 8;

    const int MAX_BLOCK_SIZE = 64 * 1024;

    // Make sure there is enough space in the buffer for compressed and uncompressed data.
    resize(stream._compressedBlock, MAX_BLOCK_SIZE);
    resize(stream._uncompressedBlock, MAX_BLOCK_SIZE);

    char * buffer = &stream._compressedBlock[0];
    int bufferSize = length(stream._compressedBlock);

    // Init gzip header
    buffer[0] = GZIP_ID1;
    buffer[1] = GZIP_ID2;
    buffer[2] = CM_DEFLATE;
    buffer[3] = FLG_FEXTRA;
    buffer[4] = 0;  // mtime
    buffer[5] = 0;
    buffer[6] = 0;
    buffer[7] = 0;
    buffer[8] = 0;
    buffer[9] = OS_UNKNOWN;
    buffer[10] = BGZF_XLEN;
    buffer[11] = 0;
    buffer[12] = BGZF_ID1;
    buffer[13] = BGZF_ID2;
    buffer[14] = BGZF_LEN;
    buffer[15] = 0;
    buffer[16] = 0;  // Placeholder for block length.
    buffer[17] = 0;

    // Loop to retry for blocks that do not compress enough.
    int inputLength = blockLength;
    int compressedLength = 0;
    while (true)
    {
        z_stream zs;
        zs.zalloc = NULL;
        zs.zfree = NULL;
        zs.next_in = static_cast<Bytef *>(static_cast<void *>(&stream._uncompressedBlock[0]));
        zs.avail_in = inputLength;
        zs.next_out = static_cast<Bytef *>(static_cast<void *>(&buffer[BLOCK_HEADER_LENGTH]));
        zs.avail_out = bufferSize - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;

        int status = deflateInit2(&zs, stream._compressLevel, Z_DEFLATED,
                                  GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY);
        if (status != Z_OK)
            return -1;  // deflateInit2() failed.

        status = deflate(&zs, Z_FINISH);
        if (status != Z_STREAM_END)
        {
            status = deflateEnd(&zs);
            if (status == Z_OK)
            {
                // Not enough space in buffer.  This can happen in the rare case the input doesn't compress enough.  We
                // try to resolve this by reducing the amount of input until it fits.
                inputLength -= 1024;
                if (inputLength <= 0)
                    return -1;  // Input reduction failed.
                continue;
            }
            return -1;  // Deflate failed.
        }

        status = deflateEnd(&zs);
        if (status != Z_OK)
            return -1;  // Deflate end failed.

        compressedLength = zs.total_out;
        compressedLength += BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
        if (compressedLength > MAX_BLOCK_SIZE)
            return -1;  // Deflate overflow, should never happen.

        break;
    }

    // Set compressed length into buffer, compute CRC and write CRC into buffer.
    _bgzfPackInt16((unsigned char*)&buffer[16], compressedLength - 1);
    __uint32 crc = crc32(0L, NULL, 0L);
    crc = crc32(crc, static_cast<Bytef *>(static_cast<void *>(&stream._uncompressedBlock[0])), inputLength);
    _bgzfPackInt32((unsigned char*)&buffer[compressedLength - 8], crc);
    _bgzfPackInt32((unsigned char*)&buffer[compressedLength - 4], inputLength);

    // Copy data that did not fit into the compressed block forward in the uncompressed data buffer.
    int remaining = blockLength - inputLength;
    if (remaining > 0)
    {
        if (remaining > inputLength)
            return -1;  // Remained too large.  Should never happen (checking here so we can use memcpy).
        memcpy(&stream._uncompressedBlock[0],
               &stream._uncompressedBlock[0] + inputLength,
               remaining);
    }
    stream._blockOffset = remaining;

    return compressedLength;
}

// ----------------------------------------------------------------------------
// Function attachToFile
// ----------------------------------------------------------------------------

// TODO(holtgrew): Rename to "reopen"? Is not not a better term, I guess.

/**
.Function.attachToFile
..class:Spec.BGZF Stream
..cat:Input/Output
..summary:Attach to already open input / output file.
..signature:attachToFile(stream, fileHandle, mode)
..param.stream:Stream to attach to file.
...type:Spec.BGZF Stream
..param.fileHandle:The file handle to attach to.
...type:nolink:$int$
..param.mode:The mode flag the file was opened with.
...type:nolink:$int$
..remarks:A previously opened file is closed.
..include:seqan/stream.
*/

inline void
attachToFile(Stream<Bgzf> & stream, int fileHandle, int mode)
{
    close(stream);  // Close previously opened file.
    stream._file.handle = fileHandle;
    stream._openMode = mode;
    stream._fileOwned = false;
    // Set compression level.
    if (mode & OPEN_WRONLY)
        stream._compressLevel = Z_DEFAULT_COMPRESSION;
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

/**
.Function.open
..class:Spec.BGZF Stream
..param.stream.type:Spec.BGZF Stream
..param.mode.remarks:When opening $Stream<Bgzf>$ for writing, you can append $'0'$-$'9'$ or $'u'$ to $"w"$ for setting the compression level to 0-9 or "no compression".
 */

inline bool
open(Stream<Bgzf> & stream, char const * filename, char const * mode)
{
    // Reset stream state.
    stream._error = 0;
    stream._openMode = 0;
    stream._blockPosition = 0;
    stream._blockLength = 0;
    stream._blockOffset = 0;
    stream._fileSize = 0;

    // Actually open files.
    if (mode[0] == 'r' || mode[0] == 'R')  // Open for reading.
    {
        stream._openMode = OPEN_RDONLY;
        open(stream._file, filename, stream._openMode);

        // Determine file size.
        if (seek(stream._file, 0, SEEK_END) == 0)
        {
            stream._error = -1;  // Seek from end of file failed.
            return false;
        }
        stream._fileSize = tell(stream._file);
        if (seek(stream._file, 0, SEEK_SET) != 0)
        {
            stream._error = -1;  // Seek from start of file failed.
            return false;
        }
    }
    else if (mode[0] == 'w' || mode[0] == 'W')  // Open for writing.
    {
        stream._compressLevel = Z_DEFAULT_COMPRESSION;

        // Set compression level.
        unsigned i;
        for (i = 0; mode[i]; ++i)
            if (mode[i] >= '0' && mode[i] <= '9') break;
        if (mode[i])
            stream._compressLevel = (int)mode[i] - '0';
        if (strchr(mode, 'u'))
            stream._compressLevel = 0;

        stream._openMode = OPEN_WRONLY | OPEN_CREATE;
        open(stream._file, filename, stream._openMode);
    }

    if (stream._file.handle != -1)  // TODO(holtgrew): File should provide isOpen().
    {
        stream._fileOwned = true;
        return true;
    }
    return false;
}

// ----------------------------------------------------------------------------
// Function streamFlush()
// ----------------------------------------------------------------------------

inline int
streamFlush(Stream<Bgzf> & stream)
{
    while (stream._blockOffset > 0)
    {
		int blockLength = _bgzfDeflateBlock(stream, stream._blockOffset);
        if (blockLength < 0)
            return -1;

        typedef Position<Stream<Bgzf> >::Type TPos;
        TPos posBefore = tell(stream._file);
        if (!write(stream._file, &stream._compressedBlock[0], blockLength))
            return -1;  // Could not write.
        TPos posAfter = tell(stream._file);
        int count = posAfter - posBefore;
        if (count != blockLength)
            return -1;  // Writing failed.

        stream._blockPosition += blockLength;
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Function close()
// ----------------------------------------------------------------------------

/**
.Function.close
..class:Spec.BGZF Stream
..signature:close(stream)
..param.stream:Stream to close.
...type:Class.Stream
 */

inline void
close(Stream<Bgzf> & stream)
{
    if (stream._file.handle == -1)
        return;   // Do nothing if file is not open.
    if (!stream._fileOwned)
        return;  // Do nothing if we do not own the file.

    if (stream._openMode & OPEN_WRONLY)
    {
        // Flush stream, compressing and writing out any uncompressed data.
        if (streamFlush(stream) != 0)
        {
            close(stream._file);
            return;  // Could not flush.
        }

        // Write an empty block.
        int blockLength = _bgzfDeflateBlock(stream, 0);
        typedef Position<Stream<Bgzf> >::Type TPos;
        TPos posBefore = tell(stream._file);
        if (!write(stream._file, &stream._compressedBlock[0], blockLength))
            return;  // Could not write
        TPos posAfter = tell(stream._file);
        int count = posAfter - posBefore;
        // if (count != blockLength)
        //     SEQAN_FAIL("Failed to write empty block at end of BGZF file.");
        (void)count;
        flush(stream._file);
    }

    // Clear the cache.
    _bgzfClearCache(stream);

    // Close file.
    close(stream._file);
}

// ----------------------------------------------------------------------------
// Function streamPeek()
// ----------------------------------------------------------------------------

inline int
streamPeek(char & c, Stream<Bgzf> & stream)
{
    // Read next block if at end of block.
    if (stream._blockOffset >= stream._blockLength)
    {
        if (_bgzfReadBlock(stream) != 0)
            return -2;  // Error.
        if (stream._blockLength == 0)
            return -1;  // EOF.
    }

    // Return next character from block without advancing stream.
    SEQAN_ASSERT_LT(stream._blockOffset, stream._blockLength);
	c = stream._uncompressedBlock[stream._blockOffset];
    return 0;
}

// ----------------------------------------------------------------------------
// Function streamReadChar()
// ----------------------------------------------------------------------------

inline int
streamReadChar(char & c, Stream<Bgzf> & stream)
{
    // Read next block if at end of block.
    if (stream._blockOffset >= stream._blockLength)
    {
        if (_bgzfReadBlock(stream) != 0)
            return -2;  // Error.
        if (stream._blockLength == 0)
            return -1;  // EOF.
    }

    // Return next address and advance stream, without reading next block yet.
	c = stream._uncompressedBlock[stream._blockOffset++];
    if (stream._blockOffset == stream._blockLength)
    {
        stream._blockPosition = tell(stream._file);
        stream._blockOffset = 0;
        stream._blockLength = 0;
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Function checkEofIsValid()
// ----------------------------------------------------------------------------

/**
.Function.checkEof
..class:Spec.BGZF Stream
..cat:Input/Output
..summary:Check that the EOF marker is present in a BGZF(/BAM) file.
..signature:streamEof(stream)
..param.stream:The BGZF Stream to check.
...type:Spec.BGZF Stream
..remarks:This is NOT equivalent to @Function.streamEof@.
..include:seqan/stream.h
*/

inline bool
checkEofIsValid(Stream<Bgzf> & stream)
{
    // Setup magic number.
	static char magic[28];
    char const * MAGIC_INIT = "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0";
    memcpy(&magic[0], MAGIC_INIT, 28);

    // Try to read magic number from end of file, then seek back to the current position.
	uint8_t buf[28];
    Position<Stream<Bgzf> >::Type offset = tell(stream._file);

    if (!seek(stream._file, -28, SEEK_END))
    {
        stream._error = -1;  // Seek from end of file failed.
        return false;
    }

    read(stream._file, &buf[0], 28);

    if (!seek(stream._file, offset, SEEK_SET))
    {
        stream._error = -1;  // Seek to original position failed.
        return false;
    }

    // Finally, compare read magic number with the expected one.
	return (memcmp(magic, buf, 28) == 0) ? true : false;
}

// ----------------------------------------------------------------------------
// Function streamEof()
// ----------------------------------------------------------------------------

inline bool
streamEof(Stream<Bgzf> & stream)
{
    if (stream._openMode & OPEN_WRONLY)
        return true;  // Always EOF when writing.
    if (stream._atEof)
        return true;  // Two consecutive streamEof() calls.

    // We are in read mode if we reach here.
    if (stream._blockOffset < stream._blockLength)
        return false;  // There is still stuff in the read buffer.
    int res =_bgzfReadBlock(stream);
    if (res == -2)
        stream._atEof = true;
    else if (res != 0)
        SEQAN_FAIL("Error reading block in streamEof().");
    else
        stream._atEof = (stream._blockLength == 0);
    return stream._atEof;
}

// TODO(holtgrew): We would rather have the const version only and then need to make some things in Stream<Bgzf> mutable.
inline bool
streamEof(Stream<Bgzf> const & stream)
{
    return streamEof(const_cast<Stream<Bgzf> &>(stream));
}

// ----------------------------------------------------------------------------
// Function streamError
// ----------------------------------------------------------------------------

inline int
streamError(Stream<Bgzf> & stream)
{
    if (stream._error == -1)
        return 0;
    return stream._error;
}

// ----------------------------------------------------------------------------
// Function streamReadBlock()
// ----------------------------------------------------------------------------

inline size_t
streamReadBlock(char * target, Stream<Bgzf> & stream, size_t maxLen)
{
    if (!(stream._openMode & OPEN_RDONLY))
        return 0;  // File not open for reading.

    // Memoize number of read bytes and pointer into the buffer target.
    size_t bytesRead = 0;
    char * destPtr = target;

    // Read at most maxLen characters, each loop iteration corresponds to reading the end of the first, the beginning of
    // the last or the whole "middle" buffers.  Of course, the first and only iteration can also only read parts of the
    // first buffer.
    while (bytesRead < maxLen)
    {
        // If there are no more bytes left in the current block then read and decompress the next block.
        int available = stream._blockLength - stream._blockOffset;
        if (available <= 0)
        {
            if (_bgzfReadBlock(stream) != 0)
                return -1;  // Could not read next block.
            available = stream._blockLength - stream._blockOffset;
            if (available <= 0)
                break;
        }

        // Copy out the number of bytes to be read or the number of available bytes in the next buffer, whichever number
        // is smaller.
        int copyLength = std::min(static_cast<int>(maxLen - bytesRead), available);
        char * buffer = &stream._uncompressedBlock[0];
        memcpy(destPtr, buffer + stream._blockOffset, copyLength);

        // Advance to next block.
        stream._blockOffset += copyLength;
        destPtr += copyLength;
        bytesRead += copyLength;
    }

    // If we read to the end of the block above then switch the block address to the next block and mark it as unread.
    if (stream._blockOffset == stream._blockLength)
    {
        stream._blockPosition = tell(stream._file);
        stream._blockOffset = 0;
        stream._blockLength = 0;
    }

    return bytesRead;
}

// ----------------------------------------------------------------------------
// Function streamWriteBlock()
// ----------------------------------------------------------------------------

inline size_t
streamWriteBlock(Stream<Bgzf> & stream, char const * source, size_t count)
{
    unsigned const MAX_BLOCK_SIZE = 64 * 1024;
    unsigned const DEFAULT_BLOCK_SIZE = 64 * 1024;

    if (!(stream._openMode & OPEN_WRONLY))
        return -1;  // File is not open for writing.

    // Allocate buffer for uncompressed data.
    resize(stream._uncompressedBlock, MAX_BLOCK_SIZE);

	char const * inPtr = source;
    int blockLength = DEFAULT_BLOCK_SIZE;
    unsigned bytesWritten = 0;

    while (bytesWritten < count)
    {
        int copyLength = std::min(static_cast<int>(blockLength - stream._blockOffset),
                                  static_cast<int>(count - bytesWritten));
        char * buffer = &stream._uncompressedBlock[0];

        memcpy(buffer + stream._blockOffset, inPtr, copyLength);

        stream._blockOffset += copyLength;
        inPtr += copyLength;
        bytesWritten += copyLength;

        if (stream._blockOffset == blockLength && streamFlush(stream) != 0)
            break;
    }

    return bytesWritten;
}

// ----------------------------------------------------------------------------
// Function streamWriteChar
// ----------------------------------------------------------------------------

inline int
streamWriteChar(Stream<Bgzf> & stream, char const c)
{
    return streamWriteBlock(stream, &c, 1) != 1;
}

// ----------------------------------------------------------------------------
// Function streamSeek()
// ----------------------------------------------------------------------------

inline int
streamSeek(Stream<Bgzf> & stream, __uint64 virtualPos, int origin)
{
    // Check whether seeking is feasible.
    if (!(stream._openMode & OPEN_RDONLY))
        return -1;  // Stream not open for reading.
    if (origin != SEEK_SET)
        return -1;  // Only seek set is implemented.

    // Extract block offset and address within the block from the virtual position.
    int blockOffset = virtualPos & 0xFFFF;
    __int64 blockAddress = (virtualPos >> 16) & 0xFFFFFFFFFFFFLL;

    // If EOF flag is set then check if we want to keep it.
    if (stream._atEof)
    {
        __int64 currentPos = tell(stream._file);
        if (currentPos != blockAddress)
            stream._atEof = false;
    }

    // Actually perform the seek.
    seek(stream._file, blockAddress, SEEK_SET);

    // Set the stream state such that the address of the block and the offset in the block are set appropriately but the
    // block is only loaded on the next read.
    stream._blockLength = 0;  // Indicates that the current block is not loaded.
    stream._blockPosition = blockAddress;
    stream._blockOffset = blockOffset;

    return 0;
}

// ----------------------------------------------------------------------------
// Function streamTell()
// ----------------------------------------------------------------------------

inline Position<Stream<Bgzf> >::Type
streamTell(Stream<Bgzf> & stream)
{
    return (stream._blockPosition << 16) | (stream._blockOffset & 0xFFFF);
}

}  // namespace seqan

#endif  // #ifndef EXTRAS_INCLUDE_SEQAN_STREAM_STREAM_BGZF_H_
