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
// A Stream specialization that works on char arrays/char *.
// ==========================================================================

// TODO(holtgrew): Using the end pointer means that the used strings are bounded which might cost some performance, although the hope is that branches are mostly predicted.

#ifndef SEQAN_STREAM_ADAPT_STREAM_CHAR_ARRAY_H_
#define SEQAN_STREAM_ADAPT_STREAM_CHAR_ARRAY_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Spec.Char Array Stream
..cat:Input/Output
..general:Class.Stream
..summary:Thin wrapper around $char *$ to the @Concept.StreamConcept|Stream@ concept.
..signature:Stream<CharArray<TPointer> >
..param.TPointer:Specification of the pointer type to work on.
...type:nolink:$char *$, $char const *$.
..remarks:This class consists of the $char *$, another $char *$ to the beginning of the array and a flag signifying EOF.
..remarks:Note that this is a bounded string variant and might have some performance problems.
..remarks:One major use case for this is to create a @Class.RecordReader@ of a string for parsing.
..include:seqan/stream.h
..example.text:Create a Char Array Stream from a @Shortcut.CharString@.
..example.code:
CharString buffer = "This is a text.";
Stream<CharArray<char const *> > stream(&buffer[0], &buffer[0] + length(buffer));

.Memfunc.Char Array Stream#Stream
..summary:Constructor
..class:Spec.Char Array Stream
..signature:CharArrayStream(ptr, ptrEnd)
..param.ptr:The $char *$ that works as the "beginning" of the stream.
..param.ptrEnd:The $char *$ that works as the "end" of the stream.
 */

template <typename TPointer>
class Stream<CharArray<TPointer> >
{
public:

    // -----------------------------------------------------------------------
    // Members
    // -----------------------------------------------------------------------
    TPointer _base;
    TPointer _ptr;
    TPointer _end;
    bool _eof;

    // -----------------------------------------------------------------------
    // Constructors
    // -----------------------------------------------------------------------

    Stream() {}

    Stream(TPointer p, TPointer e) : _base(p), _ptr(p), _end(e), _eof(false)
    {}

    Stream(Stream<CharArray<typename RemoveConst_<typename Value<TPointer>::Type>::Type *> > const & other)
            : _base(other._base), _ptr(other._ptr), _end(other._end), _eof(other._eof)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

template <>
struct Difference<Stream<CharArray<char *> > >
{
    typedef Difference<char *>::Type Type;
};

template <>
struct Difference<Stream<CharArray<char const *> > >
{
    typedef Difference<char const *>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <>
struct Position<Stream<CharArray<char *> > >
{
    typedef Position<char *>::Type Type;
};

template <>
struct Position<Stream<CharArray<char const *> > >
{
    typedef Position<char const *>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <>
struct Size<Stream<CharArray<char *> > >
{
    typedef Size<char *>::Type Type;
};

template <>
struct Size<Stream<CharArray<char const *> > >
{
    typedef Size<char const *>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <>
struct Value<Stream<CharArray<char *> > >
{
    typedef Size<char *>::Type Type;
};

template <>
struct Value<Stream<CharArray<char const *> > >
{
    typedef Value<char const *>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, IsInput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<CharArray<char *> >, IsInput>
{
    typedef True Type;
};

template <>
struct HasStreamFeature<Stream<CharArray<char const *> >, IsInput>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, IsOutput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<CharArray<char *> >, IsOutput>
{
    typedef True Type;
};

template <>
struct HasStreamFeature<Stream<CharArray<char const *> >, IsOutput>
{
    typedef False Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, HasPeek>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<CharArray<char *> >, HasPeek>
{
    typedef True Type;
};

template <>
struct HasStreamFeature<Stream<CharArray<char const *> >, HasPeek>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, HasFilename>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<CharArray<char *> >, HasFilename>
{
    typedef False Type;
};

template <>
struct HasStreamFeature<Stream<CharArray<char const *> >, HasFilename>
{
    typedef False Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, Seek<TSpec> >
// ----------------------------------------------------------------------------

template <typename TSpec>
struct HasStreamFeature<Stream<CharArray<char *> >, Seek<TSpec> >
{
    typedef True Type;
};

template <typename TSpec>
struct HasStreamFeature<Stream<CharArray<char const *> >, Seek<TSpec> >
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, Tell>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<CharArray<char *> >, Tell>
{
    typedef True Type;
};

template <>
struct HasStreamFeature<Stream<CharArray<char const *> >, Tell>
{
    typedef True Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function streamPeek()
// ----------------------------------------------------------------------------

inline int
_streamPeekImpl(char & c, char const * & ptr, char const * & end, bool & eof)
{
    if (ptr != end) {
        c = *ptr;
        return 0;
    } else {
        eof = true;
        return EOF;
    }
}

inline int
streamPeek(char & c, Stream<CharArray<char const *> > & stream)
{
    return _streamPeekImpl(c, stream._ptr, stream._end, stream._eof);
}

inline int
streamPeek(char & c, Stream<CharArray<char *> > & stream)
{
    return _streamPeekImpl(c, const_cast<char const *&>(stream._ptr), const_cast<char const *&>(stream._end), stream._eof);
}

// ----------------------------------------------------------------------------
// Function streamReadChar()
// ----------------------------------------------------------------------------

inline int
_streamReadCharImpl(bool & eof, char & c, char const * & ptr, char const * & end)
{
    if (ptr != end) {
        c = *ptr;
        ++ptr;
        return 0;
    } else {
        eof = true;
        return EOF;
    }
}

inline int
streamReadChar(char & c, Stream<CharArray<char const *> > & stream)
{
    return _streamReadCharImpl(stream._eof, c, stream._ptr, stream._end);
}

inline int
streamReadChar(char & c, Stream<CharArray<char *> > & stream)
{
    return _streamReadCharImpl(stream._eof, c, const_cast<char const *&>(stream._ptr), const_cast<char const *&>(stream._end));
}

// ----------------------------------------------------------------------------
// Function streamEof()
// ----------------------------------------------------------------------------

inline bool
streamEof(Stream<CharArray<char const *> > & stream)
{
    return stream._eof;
}

inline bool
streamEof(Stream<CharArray<char *> > & stream)
{
    return stream._eof;
}

// ----------------------------------------------------------------------------
// Function streamReadBlock()
// ----------------------------------------------------------------------------

inline size_t
_streamReadBlockImpl(bool & eof, char * target, char const * & ptr, char const * & end, size_t maxLen)
{
    char * targetPtr = target;
    size_t valuesRead = 0;
    for (size_t i = 0; i < maxLen; ++i) {
        if (ptr != end) {
            ++valuesRead;
            *targetPtr++ = *ptr++;
        } else {
            eof = true;
            break;
        }
    }
    return valuesRead;
}

inline size_t
streamReadBlock(char * target, Stream<CharArray<char const *> > & stream, size_t maxLen)
{
    return _streamReadBlockImpl(stream._eof, target, stream._ptr, stream._end, maxLen);
}

inline size_t
streamReadBlock(char * target, Stream<CharArray<char *> > & stream, size_t maxLen)
{
    return _streamReadBlockImpl(stream._eof, target, const_cast<char const *&>(stream._ptr), const_cast<char const *&>(stream._end), maxLen);
}

// ----------------------------------------------------------------------------
// Function streamWriteChar
// ----------------------------------------------------------------------------

inline int
streamWriteChar(Stream<CharArray<char *> > & stream, char const c)
{
    if (stream._ptr != stream._end) {  // Const-sized stream.
        *stream._ptr = c;
        ++stream._ptr;
        return 0;
    } else {
        stream._eof = true;
        return EOF;
    }
}

// ----------------------------------------------------------------------------
// Function streamWriteBlock()
// ----------------------------------------------------------------------------

inline size_t
streamWriteBlock(Stream<CharArray<char *> > & stream, char const * source, size_t count)
{
    size_t charsWritten = 0;
    for (; charsWritten < count && stream._ptr != stream._end; ++source, ++stream._ptr, ++charsWritten)
        *stream._ptr = *source;
    if (stream._ptr == stream._end)
        stream._eof = true;
    return charsWritten;
}

// ----------------------------------------------------------------------------
// Function streamFlush()
// ----------------------------------------------------------------------------

inline int
streamFlush(Stream<CharArray<char *> > & /*stream*/)
{
    return 0;
}

inline int
streamFlush(Stream<CharArray<char const *> > & /*stream*/)
{
    return 0;
}

// ----------------------------------------------------------------------------
// Function streamSeek()
// ----------------------------------------------------------------------------

inline int
_streamSeekImpl(char const * & ptr, char const * base, char const * end, ptrdiff_t delta, int origin)
{
    switch (origin) {
        case SEEK_SET:
            if (base + delta > end)
                return 1;
            ptr = base + delta;
            return 0;
        case SEEK_CUR:
            if (ptr + delta > end)
                return 1;
            ptr = ptr + delta;
            return 0;
        default:  // SEEK_END
            if (delta > 0)
                return 1;
            ptr = end + delta;
            return 0;
    }
}

inline int
streamSeek(Stream<CharArray<char const *> > & stream, ptrdiff_t delta, int origin)
{
    return _streamSeekImpl(stream._ptr, stream._base, stream._end, delta, origin);
}

inline int
streamSeek(Stream<CharArray<char *> > & stream, ptrdiff_t delta, int origin)
{
    return _streamSeekImpl(const_cast<char const * &>(stream._ptr), stream._base, stream._end, delta, origin);
}

// ----------------------------------------------------------------------------
// Function streamTell()
// ----------------------------------------------------------------------------

inline size_t
streamTell(Stream<CharArray<char *> > & stream)
{
    return stream._ptr - stream._base;
}

inline size_t
streamTell(Stream<CharArray<char const *> > & stream)
{
    return stream._ptr - stream._base;
}

// ----------------------------------------------------------------------------
// Function streamError
// ----------------------------------------------------------------------------

inline int
streamError(Stream<CharArray<char *> > & /*stream*/)
{
    return 0;
}

inline int
streamError(Stream<CharArray<char const *> > & /*stream*/)
{
    return 0;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_ADAPT_STREAM_CHAR_ARRAY_H_
