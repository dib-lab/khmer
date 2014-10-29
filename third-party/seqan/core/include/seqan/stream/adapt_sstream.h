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
// Adaptions for the <sstream> streams: std::istringstream,
// std::stringstream, std::ostringstream.
// ==========================================================================

#include <sstream>

#ifndef SEQAN_STREAM_ADAPT_SSTREAM_H_
#define SEQAN_STREAM_ADAPT_SSTREAM_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*
.Adaption.String stream
..cat:Input/Output
..summary:Adaption from $fstream$, $ifstream$ and $ofstream$ to the @Concept.StreamConcept@ concept.
..include:seqan/stream.h
 */

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

template <>
struct Difference< ::std::stringstream>
{
    typedef ::std::stringstream::pos_type Type;
};

template <>
struct Difference< ::std::istringstream>
{
    typedef ::std::istringstream::pos_type Type;
};

template <>
struct Difference< ::std::ostringstream>
{
    typedef ::std::ostringstream::pos_type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <>
struct Size< ::std::stringstream>
{
    typedef ::std::stringstream::pos_type Type;
};

template <>
struct Size< ::std::istringstream>
{
    typedef ::std::istringstream::pos_type Type;
};

template <>
struct Size< ::std::ostringstream>
{
    typedef ::std::ostringstream::pos_type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <>
struct Value< ::std::stringstream>
{
    typedef ::std::stringstream::char_type Type;
};

template <>
struct Value< ::std::istringstream>
{
    typedef ::std::istringstream::char_type Type;
};

template <>
struct Value< ::std::ostringstream>
{
    typedef ::std::ostringstream::char_type Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, IsInput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature< ::std::stringstream, IsInput>
{
    typedef True Type;
};

template <>
struct HasStreamFeature< ::std::istringstream, IsInput>
{
    typedef True Type;
};

template <>
struct HasStreamFeature< ::std::ostringstream, IsInput>
{
    typedef False Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature< ::std::, IsOutput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature< ::std::stringstream, IsOutput>
{
    typedef True Type;
};

template <>
struct HasStreamFeature< ::std::istringstream, IsOutput>
{
    typedef False Type;
};

template <>
struct HasStreamFeature< ::std::ostringstream, IsOutput>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature< ::std::, HasPeek>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature< ::std::stringstream, HasPeek>
{
    typedef True Type;
};

template <>
struct HasStreamFeature< ::std::istringstream, HasPeek>
{
    typedef True Type;
};

template <>
struct HasStreamFeature< ::std::ostringstream, HasPeek>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature< ::std::, HasFilename>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature< ::std::stringstream, HasFilename>
{
    typedef False Type;
};

template <>
struct HasStreamFeature< ::std::istringstream, HasFilename>
{
    typedef False Type;
};

template <>
struct HasStreamFeature< ::std::ostringstream, HasFilename>
{
    typedef False Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature< ::std::, Seek<TSpec> >
// ----------------------------------------------------------------------------

template <typename TSpec>
struct HasStreamFeature< ::std::stringstream, Seek<TSpec> >
{
    typedef True Type;
};

template <typename TSpec>
struct HasStreamFeature< ::std::istringstream, Seek<TSpec> >
{
    typedef True Type;
};

template <typename TSpec>
struct HasStreamFeature< ::std::ostringstream, Seek<TSpec> >
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature< ::std::, Tell>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature< ::std::stringstream, Tell>
{
    typedef True Type;
};

template <>
struct HasStreamFeature< ::std::istringstream, Tell>
{
    typedef True Type;
};

template <>
struct HasStreamFeature< ::std::ostringstream, Tell>
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
streamPeek(char & c, ::std::stringstream & stream)
{
    return _streamPeekImplIOStream(c, stream);
}

inline int
streamPeek(char & c, ::std::istringstream & stream)
{
    return _streamPeekImplIOStream(c, stream);
}

// ----------------------------------------------------------------------------
// Function streamReadChar()
// ----------------------------------------------------------------------------

inline int
streamReadChar(char & c, ::std::stringstream & stream)
{
    int res = _streamReadCharImplIOStream(c, stream);
    stream.seekp(stream.tellg());
    return res;
}

inline int
streamReadChar(char & c, ::std::istringstream & stream)
{
    return _streamReadCharImplIOStream(c, stream);
}

// ----------------------------------------------------------------------------
// Function streamEof()
// ----------------------------------------------------------------------------

inline bool
streamEof(::std::stringstream & stream)
{
    return stream.eof();
}

inline bool
streamEof(::std::istringstream & stream)
{
    return stream.eof();
}

inline bool
streamEof(::std::ostringstream & stream)
{
    return stream.eof();
}

// ----------------------------------------------------------------------------
// Function streamError
// ----------------------------------------------------------------------------

inline int
streamError(::std::stringstream & stream)
{
    // If we read beyond the last char, then eof and fail will be
    // set. Workaround to get the logic right.
    if ((stream.rdstate() & ::std::ios_base::eofbit) != 0)
        return 0;
    return stream.rdstate();
}

inline int
streamError(::std::istringstream & stream)
{
    // If we read beyond the last char, then eof and fail will be
    // set. Workaround to get the logic right.
    if ((stream.rdstate() & ::std::ios_base::eofbit) != 0)
        return 0;
    return stream.rdstate();
}

inline int
streamError(::std::ostringstream & stream)
{
    // If we read beyond the last char, then eof and fail will be
    // set. Workaround to get the logic right.
    if ((stream.rdstate() & ::std::ios_base::eofbit) != 0)
        return 0;
    return stream.rdstate();
}

// ----------------------------------------------------------------------------
// Function streamReadBlock()
// ----------------------------------------------------------------------------

inline size_t
streamReadBlock(char * target, ::std::stringstream & stream, size_t maxLen)
{
    int res = _streamReadBlockImplIOStream(target, stream, maxLen);
    stream.seekp(stream.tellg());
    return res;
}

inline size_t
streamReadBlock(char * target, ::std::istringstream & stream, size_t maxLen)
{
    return _streamReadBlockImplIOStream(target, stream, maxLen);
}

// ----------------------------------------------------------------------------
// Function streamWriteChar
// ----------------------------------------------------------------------------

inline int
streamWriteChar(::std::stringstream & stream, char const c)
{
    int res = _streamWriteCharImplIOStream(stream, c);
    stream.seekg(stream.tellp());
    return res;
}

inline int
streamWriteChar(::std::ostringstream & stream, char const c)
{
    return _streamWriteCharImplIOStream(stream, c);
}

// ----------------------------------------------------------------------------
// Function streamWriteBlock()
// ----------------------------------------------------------------------------

inline size_t
streamWriteBlock(::std::stringstream & stream, char const * source, size_t count)
{
    int res = _streamWriteBlockImplIOStream(stream, source, count);
    stream.seekg(stream.tellp());
    return res;
}

inline size_t
streamWriteBlock(::std::ostringstream & stream, char const * source, size_t count)
{
    return _streamWriteBlockImplIOStream(stream, source, count);
}

// ----------------------------------------------------------------------------
// Function streamPut()
// ----------------------------------------------------------------------------

// sstream
// --- characters

inline int
streamPut(::std::stringstream & stream, char const c)
{
    return streamWriteChar(stream, c);
}

template <typename TValue, typename TSpec>
inline int
streamPut(::std::stringstream & stream,
          SimpleType<TValue, TSpec> const & c)
{
    return streamWriteChar(stream, convert<char>(c));
}

// --- strings

inline int
streamPut(::std::stringstream & stream, char const * source)
{
    return (streamWriteBlock(stream, source, strlen(source))
                == strlen(source) )  ?   0 : 1;
}

template <typename TSpec>
inline int
streamPut(::std::stringstream & stream, String<char, TSpec> const & source)
{
    return (streamWriteBlock(stream, toCString(source), length(source))
                == length(source))  ?   0 : 1;
}

template <typename TValue, typename TSpec, typename TSpec2>
inline int
streamPut(::std::stringstream & stream,
          String<SimpleType<TValue, TSpec>, TSpec2> const & source)
{
    String<char, CStyle> buf = source;
    return (streamWriteBlock(stream, toCString(buf), length(buf))
                == length(buf))  ?   0 : 1;
}

// --- wildcard

template <typename TSource>
inline int
streamPut(::std::stringstream & stream, TSource const & source)
{
    stream << source;
    return stream.fail();
}

// osstream
// --- characters

inline int
streamPut(::std::ostringstream & stream, char const c)
{
    return streamWriteChar(stream, c);
}

template <typename TValue, typename TSpec>
inline int
streamPut(::std::ostringstream & stream,
          SimpleType<TValue, TSpec> const & c)
{
    return streamWriteChar(stream, convert<char>(c));
}

// --- strings

inline int
streamPut(::std::ostringstream & stream, char const * source)
{
    return (streamWriteBlock(stream, source, strlen(source))
                == strlen(source) )  ?   0 : 1;
}

template <typename TSpec>
inline int
streamPut(::std::ostringstream & stream, String<char, TSpec> const & source)
{
    return (streamWriteBlock(stream, toCString(source), length(source))
                == length(source))  ?   0 : 1;
}

// --- wildcard

template <typename TSource>
inline int
streamPut(::std::ostringstream & stream, TSource const & source)
{
    stream << source;
    return stream.fail();
}

// ----------------------------------------------------------------------------
// Function streamFlush()
// ----------------------------------------------------------------------------

inline int
streamFlush(::std::stringstream & stream)
{
    stream << std::flush;
    return 0;
}

inline int
streamFlush(::std::ostringstream & stream)
{
    stream << std::flush;
    return 0;
}

// ----------------------------------------------------------------------------
// Function streamSeek()
// ----------------------------------------------------------------------------

inline int
streamSeek(::std::stringstream & stream, long int delta, int origin)
{
    // For sstream, the input and output pointer are kept in sync.
    if (origin == SEEK_SET) {
        stream.seekg(delta, ::std::stringstream::beg);
        stream.seekp(delta, ::std::stringstream::beg);
    } else if (origin == SEEK_CUR) {
        stream.seekg(delta, ::std::stringstream::cur);
        stream.seekp(delta, ::std::stringstream::cur);
    } else {
        stream.seekg(delta, ::std::stringstream::end);
        stream.seekp(delta, ::std::stringstream::end);
    }
    return (stream.fail() || stream.bad());
}

inline int
streamSeek(::std::istringstream & stream, long int delta, int origin)
{
    if (origin == SEEK_SET)
        stream.seekg(delta, ::std::stringstream::beg);
    else if (origin == SEEK_CUR)
        stream.seekg(delta, ::std::stringstream::cur);
    else
        stream.seekg(delta, ::std::stringstream::end);
    return (stream.fail() || stream.bad());
}

inline int
streamSeek(::std::ostringstream & stream, long int delta, int origin)
{
    if (origin == SEEK_SET)
        stream.seekp(delta, ::std::stringstream::beg);
    else if (origin == SEEK_CUR)
        stream.seekp(delta, ::std::stringstream::cur);
    else
        stream.seekp(delta, ::std::stringstream::end);
    return (stream.fail() || stream.bad());
}

// ----------------------------------------------------------------------------
// Function streamTell()
// ----------------------------------------------------------------------------

inline Position< ::std::stringstream>::Type
streamTell(::std::stringstream & stream)
{
    std::streampos x = stream.tellp();
    return x;
}

inline Position< ::std::stringstream>::Type
streamTell(::std::istringstream & stream)
{
    std::streampos x = stream.tellg();
    return x;
}

inline Position< ::std::stringstream>::Type
streamTell(::std::ostringstream & stream)
{
    std::streampos x = stream.tellp();
    return x;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_ADAPT_SSTREAM_H_
