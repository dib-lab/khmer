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
// Wrapper for zlib streams.  We use a wrapper because gzFile is an alias
// to void* which could be used in other places, too.
// ==========================================================================

#include <zlib.h>

#ifndef SEQAN_STREAM_STREAM_GZ_FILE_H_
#define SEQAN_STREAM_STREAM_GZ_FILE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <> class Stream<GZFile>;
inline void close(Stream<GZFile> & stream);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Spec.GZ File Stream
..cat:Input/Output
..signature:Stream<GZFile>
..general:Class.Stream
..summary:Adaption from $gzFile$ of $<zlib.h>$ to streams.
..remarks:This is only available if @Macro.SEQAN_HAS_ZLIB@ is set to 1.
..remarks:Not copy constructable, not assignable.
..remarks:Follows the RIIA pattern when the file is opened through @Function.open@ (and thus the underlying $gzFile$ is owned).
..remarks:
Can be used as a wrapper around a $gzFile$ or create such an object itself through @Function.open@.
Also see @Memfunc.GZ File Stream#Stream|the constructor@.
..include:seqan/stream.h
..example.text:It is easy to open a GZ file via @Function.open@.
..example.code:
#include <seqan/stream.h>

Stream<GZFile> gzStream;
open(gzStream, "/path/to/file.txt.gz", "r");  // binary is implicit

// Now, work with gzStream.  The object will close the file on destruction.
..example.text:
You can also use GZ File Stream as a wrapper around an $gzFile$ object.
In this case, we have to deal with the verbose code for opening the $gzFile$ object.
..example.code:
#include <zlib.h>
#include <seqan/stream.h>

gzFile gzOut = gzopen("/path/to/file.txt.gz", "wb");
SEQAN_ASSERT(gzOut != NULL);
Stream<GZFile> gzStream(gzOut);

// Now, you can work with the stream gzStream.

// Note that you only have to close gzOut and not gzStream.
int res = gzclose(gzOut);
SEQAN_ASSERT_EQ(res, 0);

.Memfunc.GZ File Stream#Stream
..summary:Constructor
..signature:Stream()
..signature:Stream(gzFile)
..param.gzFile:The $gzFile$ to wrap.
...type:nolink:$gzFile$ (from zlib.h).
..remarks:When $gzFile$ is given then the GZ File Stream object does not own the underlying $gzFile$ and will serve as a simple wrapper.
 */

template <>
class Stream<GZFile>
{
public:
    bool _gzFileOwned;
    gzFile _gzFile;

    Stream() : _gzFileOwned(false), _gzFile(0) {}

    Stream(gzFile & gzFile) : _gzFileOwned(false), _gzFile(gzFile) {}

    ~Stream()
    {
        if (this->_gzFileOwned)
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
struct Difference<Stream<GZFile> >
{
    typedef z_off_t Type;
};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <>
struct Position<Stream<GZFile> >
{
    typedef z_off_t Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <>
struct Size<Stream<GZFile> >
{
    typedef z_off_t Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <>
struct Value<Stream<GZFile> >
{
    typedef char Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, IsInput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<GZFile>, IsInput>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, IsOutput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<GZFile>, IsOutput>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, HasPeek>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<GZFile>, HasPeek>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, HasFilename>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<GZFile>, HasFilename>
{
    typedef False Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, Seek<TSpec> >
// ----------------------------------------------------------------------------

template <typename TSpec>
struct HasStreamFeature<Stream<GZFile>, Seek<TSpec> >
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, Tell>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<GZFile>, Tell>
{
    typedef True Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

/**
.Function.open
..class:Class.Stream
..signature:open(stream, fileName, mode)
..param.stream:Stream to open.
...type:Class.Stream
..param.mode:Mode string for opening the file, e.g. $"r"$, $"w"$.
...remarks:Do not specify for @Spec.Char Array Stream@.
...type:nolink:$char const *$
 */

inline bool
open(Stream<GZFile> & stream, char const * filename, char const * mode)
{
    if (stream._gzFileOwned)
        close(stream);
    // TODO(holtgrew): Use constants instead of 0/1 for stdin/stdout.  A bit tricky to do such that it can be ported to Windows.
    if (CharString(filename) == "-")  // stdin/stdout
    {
        int fid = 0;  // stdin
        for (char const * ptr = mode; *ptr != '\0'; ++ptr)
            if (*ptr == 'w')
                fid = 1;  // stdout
        stream._gzFile = gzdopen(fid, mode);  // attach to stdout/stdin.
    }
    else
    {
        stream._gzFile = gzopen(filename, mode);
        stream._gzFileOwned = true;
    }
    if (stream._gzFile == 0)
    {
        stream._gzFileOwned = false;
        return false;
    }
    return true;
}

// ----------------------------------------------------------------------------
// Function isDirect()
// ----------------------------------------------------------------------------

/**
.Function.isDirect
..class:Spec.GZ File Stream
..cat:Input/Output
..summary:Query a GZ File Stream for being "direct."
..signature:isDirect(gzStream)
..param.gzStream:GZ File Stream to query.
...type:Spec.GZ File Stream
..returns:$bool$, indicating whether the file is opened uncompressed ("direct").
..remarks:
GZ File Streams can also open uncompressed files (at a possible performance cost).
This function returns whether the underlying file was not a compressed file and thus the file is read directly.
..include:seqan/stream.h
 */

inline bool
isDirect(Stream<GZFile> & stream)
{
    return gzdirect(stream._gzFile);
}

// ----------------------------------------------------------------------------
// Function close()
// ----------------------------------------------------------------------------

/**
.Function.close
..class:Class.Stream
..signature:close(stream)
..param.stream:Stream to close.
...type:Class.Stream
 */

inline void
close(Stream<GZFile> & stream)
{
    if (stream._gzFile == 0)
        return;
    int res = gzclose(stream._gzFile);
    (void) res;
    stream._gzFile = 0;
    return;
}

// ----------------------------------------------------------------------------
// Function streamPeek()
// ----------------------------------------------------------------------------

inline int
streamPeek(char & c, Stream<GZFile> & stream)
{
    int x = gzgetc(stream._gzFile);
    if (x < 0)
        return x;
    c = x;
    gzungetc(c, stream._gzFile);
    return 0;
}

// ----------------------------------------------------------------------------
// Function streamReadChar()
// ----------------------------------------------------------------------------

inline int
streamReadChar(char & c, Stream<GZFile> & stream)
{
    int x = gzgetc(stream._gzFile);
    if (x < 0)
        return x;
    c = x;
    return 0;
}

// ----------------------------------------------------------------------------
// Function streamEof()
// ----------------------------------------------------------------------------

inline bool
streamEof(Stream<GZFile> & stream)
{
    // std::cerr << "gzeof(stream._gzFile) == " << gzeof(stream._gzFile) << std::endl;
    return gzeof(stream._gzFile) != 0;
}

// ----------------------------------------------------------------------------
// Function streamError
// ----------------------------------------------------------------------------

inline int
streamError(Stream<GZFile> & stream)
{
    int result = 0;
    gzerror(stream._gzFile, &result);
    // Anything >= is OK/stream end etc.
    return result < 0;
}

// ----------------------------------------------------------------------------
// Function streamReadBlock()
// ----------------------------------------------------------------------------

inline size_t
streamReadBlock(char * target, Stream<GZFile> & stream, size_t maxLen)
{
    int res = gzread(stream._gzFile, target, maxLen);
    if (res < 0)
        return 0;
    return res;
}

// ----------------------------------------------------------------------------
// Function streamWriteChar
// ----------------------------------------------------------------------------

inline int
streamWriteChar(Stream<GZFile> & stream, char const c)
{
    return gzwrite(stream._gzFile, &c, 1) != 1;
}

// ----------------------------------------------------------------------------
// Function streamWriteBlock()
// ----------------------------------------------------------------------------

inline size_t
streamWriteBlock(Stream<GZFile> & stream, char const * source, size_t count)
{
    return gzwrite(stream._gzFile, const_cast<char const *>(source), count);
}

// ----------------------------------------------------------------------------
// Function streamFlush()
// ----------------------------------------------------------------------------

///.Function.streamFlush.remarks:If the stream is of type @Spec.GZ File Stream@ then this function calls $gzflush()$ with $Z_SYNC_FLUSH$. Note that many flush calls to such compressed streams reduce the compression rate.

inline int
streamFlush(Stream<GZFile> & stream)
{
    return gzflush(stream._gzFile, Z_SYNC_FLUSH);
}

// ----------------------------------------------------------------------------
// Function streamSeek()
// ----------------------------------------------------------------------------

inline int
streamSeek(Stream<GZFile> & stream, long int delta, int origin)
{
    SEQAN_ASSERT_NEQ(origin, SEEK_END);  // Not supported.
    return gzseek(stream._gzFile, delta, origin) < 0;
}

// ----------------------------------------------------------------------------
// Function streamTell()
// ----------------------------------------------------------------------------

inline Position<Stream<GZFile> >::Type
streamTell(Stream<GZFile> & stream)
{
    return gztell(stream._gzFile);
}

}  // namespace seqean

#endif  // #ifndef SEQAN_STREAM_STREAM_GZ_FILE_H_
