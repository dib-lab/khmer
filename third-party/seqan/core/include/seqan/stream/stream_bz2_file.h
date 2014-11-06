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
// Wrapper to BZFILE * that fulfills the Stream concept.  Note that, different
// from zlib, bzlib does not support seek, tell or position on BZFILE*.  We do
// not emulate support for this either.  The main use case in SeqAn is
// sequential reading and writing and we have full support for this.  We do
// not support peek either, however.  There is no flush(), it is a null
// implementation but should never be required for the use cases of bzlib in
// SeqAn.
// ==========================================================================

#include <bzlib.h>

#ifndef SEQAN_STREAM_BZ2_FILE_WRAPPER_H_
#define SEQAN_STREAM_BZ2_FILE_WRAPPER_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <> class Stream<BZ2File>;
inline void close(Stream<BZ2File> & stream);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Spec.BZ2 File Stream
..cat:Input/Output
..signature:Stream<BZ2File>
..general:Class.Stream
..summary:Wrapper for $BZFILE *$ streams from bzlib.
..remarks:This is only available if @Macro.SEQAN_HAS_ZLIB@ is set to 1.
..remarks:Not copy constructable.
..remarks:Follows the RIIA pattern when the file is opened through @Function.open@ (and thus the underlying $BZFILE *$ and $FILE *$ are owned).
..remarks:
Can be used as a wrapper around a $BZFILE *$ or create such an object itself through @Function.open@.
Also see @Memfunc.BZ2 File Stream#Stream|the constructor@.
..include:seqan/stream.h
..example.text:It is easy to open a BZ2 file via @Function.open@.
..example.code:
#include <seqan/stream.h>

Stream<BZ2File> bzStream;
open(bzStream, "/path/to/file.txt.bz2", "r");  // binary is implicit

// Now, work with bzStream.  The object will close the file on destruction.
..example.text:
You can also use BZ2 File Stream as a wrapper around $BZFILE *$.
In this case, we have to deal with the verbose code for opening $FILE *$ and $BZFILE *$ by hand.
..example.code:
#include <cstdio>
#include <bzlib.h>
#include <seqan/stream.h>

// Open normal FILE *, BZFILE * above, and Stream<BZ2File> on top.
FILE * f = fopen("/path/to/file.txt.bz2", "rb");
SEQAN_ASSERT(f != NULL);
int err = BZ_OK;
BZFILE * f2 = BZ2_bzReadOpen(&err, f, 0, 0, NULL, 0);
SEQAN_ASSERT_EQ(err, BZ_OK);
Stream<BZ2File> bzStream(f2);

// Now, you can work with the stream bzStream.

// Note that you only have to close f2 and f, not bzStream.
BZ2_bzReadClose(&err, f2);
SEQAN_ASSERT_EQ(err, BZ_OK);
fclose(f);

.Memfunc.BZ2 File Stream#Stream
..summary:Constructor
..class:Spec.BZ2 File Stream
..signature:Stream<BZ2File>()
..signature:Stream<BZ2File>(file)
..param.file:The $BZFILE *$ to wrap.
...type:nolink:$BZFILE *$.
..remarks:When $file$ is given then the BZ2 File Stream object does not own the underlying $BZFILE *$ object and will serve as a simple wrapper.
 */

template <>
class Stream<BZ2File>
{
public:
    bool _fileOwned;
    BZFILE * _file;
    FILE * _underlyingFile;
    int _error;
    char _rw;

    Stream() : _fileOwned(false), _file(0), _underlyingFile(0), _error(0), _rw('-')
    {}

    Stream(BZFILE * file) : _fileOwned(false), _file(file), _error(0), _rw('-')
    {}

    ~Stream()
    {
        if (this->_fileOwned)
            close(*this);
    }

private:
    // Disable default, copy construction and assignment.
    Stream(Stream const & /*other*/) {}
    Stream & operator=(Stream const & /*other*/) { return *this; }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

template <>
struct Difference<Stream<BZ2File> >
{
    // bzlib streams rely on FILE * streams
    typedef Difference<Stream<FILE *> >::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <>
struct Position<Stream<BZ2File> >
{
    // bzlib streams rely on FILE * streams
    typedef Position<Stream<FILE *> >::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <>
struct Size<Stream<BZ2File> >
{
    // bzlib streams rely on FILE * streams
    typedef Size<Stream<FILE *> >::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <>
struct Value<Stream<BZ2File> >
{
    // bzlib streams rely on FILE * streams
    typedef Value<Stream<FILE *> >::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, IsInput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<BZ2File>, IsInput>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, IsOutput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<BZ2File>, IsOutput>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, HasPeek>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<BZ2File>, HasPeek>
{
    typedef False Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, HasFilename>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<BZ2File>, HasFilename>
{
    typedef False Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, Seek<TSpec> >
// ----------------------------------------------------------------------------

template <typename TSpec>
struct HasStreamFeature<Stream<BZ2File>, Seek<TSpec> >
{
    typedef False Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, Tell>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<BZ2File>, Tell>
{
    typedef False Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

inline bool
open(Stream<BZ2File> & stream, char const * filename, char const * mode)
{
    if (stream._fileOwned)
        close(stream);
    CharString modeStr = mode;
    if (length(modeStr) == 0u || (modeStr[0] != 'r' && modeStr[0] != 'w'))
        return false;
    if (modeStr == "r" || modeStr == "w")
        appendValue(modeStr, 'b');
    stream._rw = modeStr[0];
    if (CharString(filename) == "-")
    {
        stream._fileOwned = false;
        // TODO(holtgrew): Use constants instead of 0/1 for stdin/stdout.  A bit tricky to do such that it can be ported to Windows.
        if (stream._rw == 'r')
            stream._underlyingFile = stdin;
        else
            stream._underlyingFile = stdout;
    }
    else
    {
        stream._underlyingFile = fopen(filename, toCString(modeStr));
        if (stream._underlyingFile == 0)
            return false;
        stream._fileOwned = true;
    }
    if (stream._rw == 'w')
        stream._file = BZ2_bzWriteOpen(&stream._error, stream._underlyingFile, 7, 0, 0);
    else
        stream._file = BZ2_bzReadOpen(&stream._error, stream._underlyingFile, 0, 0, NULL, 0);
    if (stream._file == 0 || stream._error != 0)
    {
        stream._file = 0;
        stream._underlyingFile = 0;
        stream._fileOwned = false;
        return false;
    }
    return true;
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
close(Stream<BZ2File> & stream)
{
    if (stream._file == 0 || stream._file == 0)
        return;
    if (stream._rw == 'w')
        BZ2_bzWriteClose(&stream._error, stream._file, 0, NULL, NULL);
    else
        BZ2_bzReadClose(&stream._error, stream._file);
    fclose(stream._underlyingFile);
    stream._file = 0;
    stream._underlyingFile = 0;
    stream._rw = '-';
    stream._fileOwned = false;
}

// ----------------------------------------------------------------------------
// Function streamEof()
// ----------------------------------------------------------------------------

inline bool
streamEof(Stream<BZ2File> & stream)
{
    // std::cerr << "stream._error == " << stream._error << std::endl;
    // std::cerr << "BZ_FINISH == " << BZ_STREAM_END << std::endl;
    return stream._error == BZ_STREAM_END;
}

// ----------------------------------------------------------------------------
// Function streamTell()
// ----------------------------------------------------------------------------

// Always returns 0.

inline Position<Stream<BZ2File> >::Type
streamTell(Stream<BZ2File> const & /*stream*/)
{
    return 0;
}

// ----------------------------------------------------------------------------
// Function streamReadChar()
// ----------------------------------------------------------------------------

inline int
streamReadChar(char & c, Stream<BZ2File> & stream)
{
    if (streamEof(stream))
        return 1;
    return BZ2_bzRead(&stream._error, stream._file, &c, 1) != 1;
}

// ----------------------------------------------------------------------------
// Function streamError
// ----------------------------------------------------------------------------

inline int
streamError(Stream<BZ2File> & stream)
{
    // Anything >= means OK.
    if (stream._error < 0)
        return stream._error;
    return 0;
}

// ----------------------------------------------------------------------------
// Function streamReadBlock()
// ----------------------------------------------------------------------------

inline size_t
streamReadBlock(char * target, Stream<BZ2File> & stream, size_t maxLen)
{
    return BZ2_bzRead(&stream._error, stream._file, target, maxLen);
}

// ----------------------------------------------------------------------------
// Function streamWriteChar
// ----------------------------------------------------------------------------

inline int
streamWriteChar(Stream<BZ2File> & stream, char const c)
{
    BZ2_bzWrite(&stream._error, stream._file, const_cast<char *>(&c), 1);
    return stream._error;
}

// ----------------------------------------------------------------------------
// Function streamWriteBlock()
// ----------------------------------------------------------------------------

inline size_t
streamWriteBlock(Stream<BZ2File> & stream, char const * source, size_t count)
{
    BZ2_bzWrite(&stream._error, stream._file, const_cast<char *>(source), count);
    if (stream._error)
        return 0;
    else
        return count;
}

// ----------------------------------------------------------------------------
// Function streamFlush()
// ----------------------------------------------------------------------------

///.Function.streamFlush.remarks:If the stream is of type @Spec.BZ2 File Stream@ then this function does nothing.

inline int
streamFlush(Stream<BZ2File> & /*stream*/)
{
    // Null implementation, function is not so important in the use cases for
    // SeqAn anyway.
    return 0;
}

}  // namespace seqean

#endif  // #ifndef SEQAN_STREAM_BZ2_FILE_WRAPPER_H_
