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
// Adaptions for the <fstream> streams: std::ifstream, std::fstream,
// std::ofstream.
// ==========================================================================

#include <fstream>

#ifndef SEQAN_STREAM_ADAPT_FSTREAM_H_
#define SEQAN_STREAM_ADAPT_FSTREAM_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Adaption.File stream
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
struct Difference< ::std::fstream>
{
    typedef ::std::fstream::pos_type Type;
};

template <>
struct Difference< ::std::ifstream>
{
    typedef ::std::ifstream::pos_type Type;
};

template <>
struct Difference< ::std::ofstream>
{
    typedef ::std::ofstream::pos_type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <>
struct Size< ::std::fstream>
{
    typedef size_t Type;
};

template <>
struct Size< ::std::ifstream>
{
    typedef size_t Type;
};

template <>
struct Size< ::std::ofstream>
{
    typedef size_t Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <>
struct Value< ::std::fstream>
{
    typedef ::std::fstream::char_type Type;
};

template <>
struct Value< ::std::ifstream>
{
    typedef ::std::ifstream::char_type Type;
};

template <>
struct Value< ::std::ofstream>
{
    typedef ::std::ofstream::char_type Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, IsInput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature< ::std::fstream, IsInput>
{
    typedef True Type;
};

template <>
struct HasStreamFeature< ::std::ifstream, IsInput>
{
    typedef True Type;
};

template <>
struct HasStreamFeature< ::std::ofstream, IsInput>
{
    typedef False Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature< ::std::, IsOutput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature< ::std::fstream, IsOutput>
{
    typedef True Type;
};

template <>
struct HasStreamFeature< ::std::ifstream, IsOutput>
{
    typedef False Type;
};

template <>
struct HasStreamFeature< ::std::ofstream, IsOutput>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature< ::std::, HasPeek>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature< ::std::fstream, HasPeek>
{
    typedef True Type;
};

template <>
struct HasStreamFeature< ::std::ifstream, HasPeek>
{
    typedef True Type;
};

template <>
struct HasStreamFeature< ::std::ofstream, HasPeek>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature< ::std::, HasFilename>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature< ::std::fstream, HasFilename>
{
    typedef False Type;
};

template <>
struct HasStreamFeature< ::std::ifstream, HasFilename>
{
    typedef False Type;
};

template <>
struct HasStreamFeature< ::std::ofstream, HasFilename>
{
    typedef False Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature< ::std::, Seek<TSpec> >
// ----------------------------------------------------------------------------

template <typename TSpec>
struct HasStreamFeature< ::std::fstream, Seek<TSpec> >
{
    typedef True Type;
};

template <typename TSpec>
struct HasStreamFeature< ::std::ifstream, Seek<TSpec> >
{
    typedef True Type;
};

template <typename TSpec>
struct HasStreamFeature< ::std::ofstream, Seek<TSpec> >
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature< ::std::, Tell>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature< ::std::fstream, Tell>
{
    typedef True Type;
};

template <>
struct HasStreamFeature< ::std::ifstream, Tell>
{
    typedef True Type;
};

template <>
struct HasStreamFeature< ::std::ofstream, Tell>
{
    typedef True Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function streamPeek()
// ----------------------------------------------------------------------------

template <typename TStream>
inline int
_streamPeekImplIOStream(char & c, TStream & stream)
{
    int x = stream.peek();
    if (x == EOF)
        return EOF;
    c = x;
    return 0;
}

inline int
streamPeek(char & c, ::std::fstream & stream)
{
    return _streamPeekImplIOStream(c, stream);
}

inline int
streamPeek(char & c, ::std::ifstream & stream)
{
    return _streamPeekImplIOStream(c, stream);
}

// ----------------------------------------------------------------------------
// Function streamReadChar()
// ----------------------------------------------------------------------------

template <typename TStream>
inline int
_streamReadCharImplIOStream(char & c, TStream & stream)
{
    int x = stream.get();
    if (x == EOF)
        return EOF;
    c = x;
    return 0;
}

inline int
streamReadChar(char & c, ::std::fstream & stream)
{
    return _streamReadCharImplIOStream(c, stream);
}

inline int
streamReadChar(char & c, ::std::ifstream & stream)
{
    return _streamReadCharImplIOStream(c, stream);
}

// ----------------------------------------------------------------------------
// Function streamEof()
// ----------------------------------------------------------------------------

inline bool
streamEof(::std::fstream & stream)
{
    return stream.eof();
}

inline bool
streamEof(::std::ifstream & stream)
{
    return stream.eof();
}

inline bool
streamEof(::std::ofstream & stream)
{
    return stream.eof();
}

// ----------------------------------------------------------------------------
// Function streamError
// ----------------------------------------------------------------------------

inline int
streamError(::std::fstream & stream)
{
    // If we read beyond the last char, then eof and fail will be
    // set. Workaround to get the logic right.
    if ((stream.rdstate() & ::std::ios_base::eofbit) != 0)
        return 0;
    return stream.rdstate();
}

inline int
streamError(::std::ifstream & stream)
{
    // If we read beyond the last char, then eof and fail will be
    // set. Workaround to get the logic right.
    if ((stream.rdstate() & ::std::ios_base::eofbit) != 0)
        return 0;
    return stream.rdstate();
}

inline int
streamError(::std::ofstream & stream)
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

template <typename TStream>
inline size_t
_streamReadBlockImplIOStream(char * target, TStream & stream, size_t maxLen)
{
    stream.read(target, maxLen);
    return stream.gcount();
}

inline size_t
streamReadBlock(char * target, ::std::fstream & stream, size_t maxLen)
{
    return _streamReadBlockImplIOStream(target, stream, maxLen);
}

inline size_t
streamReadBlock(char * target, ::std::ifstream & stream, size_t maxLen)
{
    return _streamReadBlockImplIOStream(target, stream, maxLen);
}

// ----------------------------------------------------------------------------
// Function streamWriteChar
// ----------------------------------------------------------------------------

template <typename TStream>
inline int
_streamWriteCharImplIOStream(TStream & stream, char const c)
{
    stream.put(c);
    if (!stream.bad())
        return 0;
    if (stream.eof())
        return EOF;
    SEQAN_ASSERT(stream.fail());
    return 1;
}

inline int
streamWriteChar(::std::fstream & stream, char const c)
{
    return _streamWriteCharImplIOStream(stream, c);
}

inline int
streamWriteChar(::std::ofstream & stream, char const c)
{
    return _streamWriteCharImplIOStream(stream, c);
}

// ----------------------------------------------------------------------------
// Function streamWriteBlock()
// ----------------------------------------------------------------------------

template <typename TStream>
inline size_t
_streamWriteBlockImplIOStream(TStream & stream, char const * source, size_t count)
{
    stream.write(source, count);
    return count; 
}

inline size_t
streamWriteBlock(::std::fstream & stream, char const * source, size_t count)
{
    return _streamWriteBlockImplIOStream(stream, source, count);
}

inline size_t
streamWriteBlock(::std::ofstream & stream, char const * source, size_t count)
{
    return _streamWriteBlockImplIOStream(stream, source, count);
}

// ----------------------------------------------------------------------------
// Function streamPut()
// ----------------------------------------------------------------------------

// fstream
// --- characters

inline int
streamPut(::std::fstream & stream, char const c)
{
    return streamWriteChar(stream, c);
}

template <typename TValue, typename TSpec>
inline int
streamPut(::std::fstream & stream,
          SimpleType<TValue, TSpec> const & c)
{
    return streamWriteChar(stream, convert<char>(c));
}

// --- strings

inline int
streamPut(::std::fstream & stream, char const * source)
{
    return (streamWriteBlock(stream, source, strlen(source))
                == strlen(source) )  ?   0 : 1;
}

template <typename TSpec>
inline int
streamPut(::std::fstream & stream, String<char, TSpec> const & source)
{
    return (streamWriteBlock(stream, toCString(source), length(source))
                == length(source))  ?   0 : 1;
}

template <typename TValue, typename TSpec, typename TSpec2>
inline int
streamPut(::std::fstream & stream,
          String<SimpleType<TValue, TSpec>, TSpec2> const & source)
{
    String<char, CStyle> buf = source;
    return (streamWriteBlock(stream, toCString(buf), length(buf))
                == length(buf))  ?   0 : 1;
}

// --- wildcard

template <typename TSource>
inline int
streamPut(::std::fstream & stream, TSource const & source)
{
    stream << source;
    return stream.fail();
}

// ofstream
// --- characters

inline int
streamPut(::std::ofstream & stream, char const c)
{
    return streamWriteChar(stream, c);
}

template <typename TValue, typename TSpec>
inline int
streamPut(::std::ofstream & stream,
          SimpleType<TValue, TSpec> const & c)
{
    return streamWriteChar(stream, convert<char>(c));
}

// --- strings

inline int
streamPut(::std::ofstream & stream, char const * source)
{
    return (streamWriteBlock(stream, source, strlen(source))
                == strlen(source) )  ?   0 : 1;
}

template <typename TSpec>
inline int
streamPut(::std::ofstream & stream, String<char, TSpec> const & source)
{
    return (streamWriteBlock(stream, toCString(source), length(source))
                == length(source))  ?   0 : 1;
}

// --- wildcard

template <typename TSource>
inline int
streamPut(::std::ofstream & stream, TSource const & source)
{
    stream << source;
    return stream.fail();
}

// ----------------------------------------------------------------------------
// Function streamFlush()
// ----------------------------------------------------------------------------

inline int
streamFlush(::std::fstream & stream)
{
    stream << std::flush;
    return 0;
}

inline int
streamFlush(::std::ofstream & stream)
{
    stream << std::flush;
    return 0;
}

// ----------------------------------------------------------------------------
// Function streamSeek()
// ----------------------------------------------------------------------------

inline int
streamSeek(::std::fstream & stream, long int delta, int origin)
{
    stream.clear();  // Reset fail flags before seek.
    // For fstream, the input and output pointer are kept in sync.
    if (origin == SEEK_SET) {
        stream.seekg(delta, ::std::fstream::beg);
    } else if (origin == SEEK_CUR) {
        stream.seekg(delta, ::std::fstream::cur);
    } else {
        stream.seekg(delta, ::std::fstream::end);
    }
    int res = (stream.fail() || stream.bad());
    stream.clear();  // Reset EOF flag.
    return res;
}

inline int
streamSeek(::std::ifstream & stream, long int delta, int origin)
{
    stream.clear();  // Reset fail flags before seek.
    if (origin == SEEK_SET)
        stream.seekg(delta, ::std::fstream::beg);
    else if (origin == SEEK_CUR)
        stream.seekg(delta, ::std::fstream::cur);
    else
        stream.seekg(delta, ::std::fstream::end);
    int res = (stream.fail() || stream.bad());
    stream.clear();  // Reset EOF flag.
    return res;
}

inline int
streamSeek(::std::ofstream & stream, long int delta, int origin)
{
    stream.clear();  // Reset fail flags before seek.
    if (origin == SEEK_SET)
        stream.seekp(delta, ::std::fstream::beg);
    else if (origin == SEEK_CUR)
        stream.seekp(delta, ::std::fstream::cur);
    else
        stream.seekp(delta, ::std::fstream::end);
    int res = (stream.fail() || stream.bad());
    stream.clear();  // Reset EOF flag.
    return res;
}

// ----------------------------------------------------------------------------
// Function streamTell()
// ----------------------------------------------------------------------------

inline Position< ::std::fstream>::Type
streamTell(::std::fstream & stream)
{
    SEQAN_ASSERT_EQ(stream.tellp(), stream.tellg());
    std::streampos x = stream.tellp();
    return x;
}

inline Position< ::std::fstream>::Type
streamTell(::std::ifstream & stream)
{
    std::streampos x = stream.tellg();
    return x;
}

inline Position< ::std::fstream>::Type
streamTell(::std::ofstream & stream)
{
    std::streampos x = stream.tellp();
    return x;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_ADAPT_FSTREAM_H_
