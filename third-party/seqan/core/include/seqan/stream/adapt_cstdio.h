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
// Adaption for <cstdio> streams: std::FILE * to the stream concept.
// ==========================================================================

#include <cstdio>
#include <sstream>

#ifndef SEQAN_STREAM_ADAPT_CSTIO_H_
#define SEQAN_STREAM_ADAPT_CSTIO_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Adaption."FILE *"
..cat:Input/Output
..summary:Adaption from $FILE *$ of $<cstdio>$ to streams.
..include:seqan/stream.h
 */

// ============================================================================
// Metafunctions
// ============================================================================

/* // Clashes with definition of these metafunctions in file module.
// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

template <>
struct Difference<FILE *>
{
    typedef long int Type;
};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <>
struct Position<FILE *>
{
    typedef long int Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <>
struct Size<FILE *>
{
    typedef long int Type;
};
*/

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, IsInput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<FILE *, IsInput>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, IsOutput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<FILE *, IsOutput>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, HasPeek>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<FILE *, HasPeek>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, HasFilename>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<FILE *, HasFilename>
{
    typedef False Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, Seek<TSpec> >
// ----------------------------------------------------------------------------

template <typename TSpec>
struct HasStreamFeature<FILE *, Seek<TSpec> >
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, Tell>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<FILE *, Tell>
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
streamPeek(char & c, FILE * stream)
{
    c = fgetc(stream);
    ungetc(c, stream);
    if (c == EOF)
        return EOF;
    return 0;
}

// ----------------------------------------------------------------------------
// Function streamReadChar()
// ----------------------------------------------------------------------------

inline int
streamReadChar(char & c, FILE * stream)
{
    c = fgetc(stream);
    if (c == EOF)
        return EOF;
    return 0;
}

// ----------------------------------------------------------------------------
// Function streamEof()
// ----------------------------------------------------------------------------

inline bool
streamEof(FILE * stream)
{
    return ::std::feof(stream) != 0;
}

// ----------------------------------------------------------------------------
// Function streamError
// ----------------------------------------------------------------------------

inline int
streamError(FILE * stream)
{
    return ::std::ferror(stream);
}

// ----------------------------------------------------------------------------
// Function streamReadBlock()
// ----------------------------------------------------------------------------

inline size_t
streamReadBlock(char * target, FILE * stream, size_t maxLen)
{
    return ::std::fread(target, sizeof(char), maxLen, stream);
}

// ----------------------------------------------------------------------------
// Function streamWriteChar
// ----------------------------------------------------------------------------

inline int
streamWriteChar(FILE * stream, char const c)
{
    int x = ::std::fputc(c, stream);
    if (x == EOF)
        return EOF;
    return c != x;
}

// ----------------------------------------------------------------------------
// Function streamWriteBlock()
// ----------------------------------------------------------------------------

inline size_t
streamWriteBlock(FILE * stream, char const * source, size_t count)
{
    return ::std::fwrite(source, sizeof(char), count, stream);
}

// ----------------------------------------------------------------------------
// Function streamPut()
// ----------------------------------------------------------------------------

// --- strings

inline int
streamPut(FILE * stream, char const * source)
{
    return (streamWriteBlock(stream, source, strlen(source))
                == strlen(source) )  ?   0 : 1;
}

template <typename TSpec>
inline int
streamPut(FILE * stream, String<char, TSpec> const & source)
{
    return (streamWriteBlock(stream, toCString(source), length(source))
                == length(source))  ?   0 : 1;
}

template <typename TSpec, typename TSpec2>
inline int
streamPut(FILE * stream,
          String<SimpleType<unsigned char, TSpec>, TSpec2> const & source)
{
    String<char, CStyle> buf = source;
    return (streamWriteBlock(stream, toCString(buf), length(buf))
                == length(buf))  ?   0 : 1;
}

// --- characters

inline int
streamPut(FILE * stream, char const c)
{
    return streamWriteChar(stream, c);
}

template <typename TValue, typename TSpec>
inline int
streamPut(FILE * stream, SimpleType<TValue, TSpec> const & c)
{
    return streamPut(stream, convert<char>(c));
}

// --- numbers

inline char const *
_streamPutChar(char const * /*tag*/)
{
    return "%s";
}

inline char const *
_streamPutChar(int const /*tag*/)
{
    return "%d";
}

inline char const *
_streamPutChar(unsigned int const /*tag*/)
{
    return "%u";
}

inline char const *
_streamPutChar(long const /*tag*/)
{
    return "%ld";
}

inline char const *
_streamPutChar(unsigned long const /*tag*/)
{
    return "%lu";
}

inline char const *
_streamPutChar(float const /*tag*/)
{
    return "%.2f"; 
}

inline char const *
_streamPutChar(double const /*tag*/)
{
    return "%.2lf";
}

// template <typename TValue, typename TSpec>
// inline char const *
// _streamPutChar(SimpleType<TValue, TSpec> const & /*tag*/)
// {
//     return _streamPutChar(TValue());
// }


// TODO(h4nn3s) according to man fprintf's point character is locale dependent,
// maybe overload for doubles and floats to avoid that?
template <typename TSource>
inline int
streamPut(FILE * stream, TSource const & source)
{
    std::stringstream tmp;
    tmp << source;
    return streamWriteBlock(stream, tmp.str().c_str(), tmp.str().size()) != tmp.str().size();
}


// ----------------------------------------------------------------------------
// Function streamFlush()
// ----------------------------------------------------------------------------

inline int
streamFlush(FILE * stream)
{
    return ::std::fflush(stream);
}

// ----------------------------------------------------------------------------
// Function streamSeek()
// ----------------------------------------------------------------------------

inline int
streamSeek(FILE * stream, long int delta, int origin)
{
    return ::std::fseek(stream, delta, origin);
}

// ----------------------------------------------------------------------------
// Function streamTell()
// ----------------------------------------------------------------------------

inline Position<FILE *>::Type
streamTell(FILE * stream)
{
    return ::std::ftell(stream);
}

}  // namespace seqean

#endif  // #ifndef SEQAN_STREAM_ADAPT_CSTIO_H_
