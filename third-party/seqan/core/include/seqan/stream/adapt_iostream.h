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
// Adaptions for <iostream> strams: std::istream, std::ostream like std::cin.
// ==========================================================================

#ifndef SEQAN_STREAM_ADAPT_IOSTREAM_H_
#define SEQAN_STREAM_ADAPT_IOSTREAM_H_

#include <iostream>

namespace seqan {

// TODO(holtgrew): Copied from adapt_fstream.h, not tested! Maybe specialize templates for basic_* streams?

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Adaption.IO stream
..summary:Adaption of standard C++ iostream objects.
..cat:Input/Output
..remarks:Also, adaption from $ostream$ and $istream$ to the @Concept.StreamConcept@ concept.
..include:seqan/stream.h
 */

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

template <>
struct Difference< ::std::ostream>
{
    typedef ::std::ostream::pos_type Type;
};

template <>
struct Difference< ::std::istream>
{
    typedef ::std::istream::pos_type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <>
struct Size< ::std::ostream>
{
    typedef ::std::ostream::pos_type Type;
};

template <>
struct Size< ::std::istream>
{
    typedef ::std::istream::pos_type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <>
struct Value< ::std::ostream>
{
    typedef ::std::ostream::char_type Type;
};

template <>
struct Value< ::std::istream>
{
    typedef ::std::istream::char_type Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, IsInput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature< ::std::ostream, IsInput>
{
    typedef False Type;
};

template <>
struct HasStreamFeature< ::std::istream, IsInput>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature< ::std::, IsOutput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature< ::std::ostream, IsOutput>
{
    typedef True Type;
};

template <>
struct HasStreamFeature< ::std::istream, IsOutput>
{
    typedef False Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature< ::std::, HasPeek>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature< ::std::ostream, HasPeek>
{
    typedef False Type;
};

template <>
struct HasStreamFeature< ::std::istream, HasPeek>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature< ::std::, HasFilename>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature< ::std::ostream, HasFilename>
{
    typedef False Type;
};

template <>
struct HasStreamFeature< ::std::istream, HasFilename>
{
    typedef False Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature< ::std::, Seek<TSpec> >
// ----------------------------------------------------------------------------

template <typename TSpec>
struct HasStreamFeature< ::std::ostream, Seek<TSpec> >
{
    typedef True Type;
};

template <typename TSpec>
struct HasStreamFeature< ::std::istream, Seek<TSpec> >
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature< ::std::, Tell>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature< ::std::ostream, Tell>
{
    typedef False Type;
};

template <>
struct HasStreamFeature< ::std::istream, Tell>
{
    typedef False Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function streamPeek()
// ----------------------------------------------------------------------------

inline int
streamPeek(char & c, ::std::istream & stream)
{
    return _streamPeekImplIOStream(c, stream);
}

// ----------------------------------------------------------------------------
// Function streamReadChar()
// ----------------------------------------------------------------------------

inline int
streamReadChar(char & c, ::std::istream & stream)
{
    return _streamReadCharImplIOStream(c, stream);
}

// ----------------------------------------------------------------------------
// Function streamEof()
// ----------------------------------------------------------------------------

inline bool
streamEof(::std::ostream & stream)
{
    return stream.eof();
}

inline bool
streamEof(::std::istream & stream)
{
    return stream.eof();
}

// ----------------------------------------------------------------------------
// Function streamError
// ----------------------------------------------------------------------------

inline int
streamError(::std::ostream & stream)
{
    // If we read beyond the last char, then eof and fail will be
    // set. Workaround to get the logic right.
    if ((stream.rdstate() & ::std::ios_base::eofbit) != 0)
        return 0;
    return stream.rdstate();
}

inline int
streamError(::std::istream & stream)
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
streamReadBlock(char * target, ::std::istream & stream, size_t maxLen)
{
    return _streamReadBlockImplIOStream(target, stream, maxLen);
}

// ----------------------------------------------------------------------------
// Function streamWriteChar
// ----------------------------------------------------------------------------

inline int
streamWriteChar(::std::ostream & stream, char const c)
{
    return _streamWriteCharImplIOStream(stream, c);
}

// ----------------------------------------------------------------------------
// Function streamWriteBlock()
// ----------------------------------------------------------------------------

inline size_t
streamWriteBlock(::std::ostream & stream, char const * source, size_t count)
{
    return _streamWriteBlockImplIOStream(stream, source, count);
}

// ----------------------------------------------------------------------------
// Function streamTell()
// ----------------------------------------------------------------------------

inline Position< ::std::istream>::Type
streamTell(::std::istream & stream)
{
    std::streampos x = stream.tellg();
    return x;
}

inline Position< ::std::ostream>::Type
streamTell(::std::ostream & stream)
{
    std::streampos x = stream.tellp();
    return x;
}

// ----------------------------------------------------------------------------
// Function streamPut()
// ----------------------------------------------------------------------------

// ostream
// --- characters

inline int
streamPut(::std::ostream & stream, char const c)
{
    return streamWriteChar(stream, c);
}

template <typename TValue, typename TSpec>
inline int
streamPut(::std::ostream & stream,
          SimpleType<TValue, TSpec> const & c)
{
    return streamWriteChar(stream, convert<char>(c));
}

// --- strings

inline int
streamPut(::std::ostream & stream, char const * source)
{
    return (streamWriteBlock(stream, source, strlen(source))
                == strlen(source) )  ?   0 : 1;
}

template <typename TSpec>
inline int
streamPut(::std::ostream & stream, String<char, TSpec> const & source)
{
    return (streamWriteBlock(stream, toCString(source), length(source))
                == length(source))  ?   0 : 1;
}

// --- wildcard

template <typename TSource>
inline int
streamPut(::std::ostream & stream, TSource const & source)
{
    stream << source;
    return stream.fail();
}

// ----------------------------------------------------------------------------
// Function streamFlush()
// ----------------------------------------------------------------------------

inline int
streamFlush(::std::ostream & stream)
{
    stream << std::flush;
    return 0;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_ADAPT_IOSTREAM_H_
