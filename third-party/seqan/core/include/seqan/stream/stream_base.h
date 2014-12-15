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
// Base class for streams.
// ==========================================================================

#ifndef SEQAN_STREAM_STREAM_BASE_H_
#define SEQAN_STREAM_STREAM_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Class.Stream
..cat:Input/Output
..signature:Stream<TSpec>
..summary:Abstract base class to fulfill the @Concept.StreamConcept@ concept.
..concept:Concept.StreamConcept
..include:seqan/stream.h
 */

template <typename TPointer = char *>
struct CharArray;

#if SEQAN_HAS_ZLIB  // Enable Stream<GZFile> if available.
struct GZFile_;
typedef Tag<GZFile_> GZFile;
#endif  // #if SEQAN_HAS_ZLIB

#if SEQAN_HAS_BZIP2  // Enable Stream<BZ2File> if available.
struct BZ2File_;
typedef Tag<BZ2File_> BZ2File;
#endif  // #if SEQAN_HAS_ZLIB

template <typename TSpec>
class Stream;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

///.Function.atEnd.param.iterator.type:Class.Stream
///.Function.atEnd.class:Class.Stream

template <typename TSpec>
inline bool
atEnd(Stream<TSpec> & stream)
{
    return streamEof(stream);
}

template <typename TSpec>
inline bool
atEnd(Stream<TSpec> const & stream)
{
    return streamEof(stream);
}

// ----------------------------------------------------------------------------
// Function streamPut()
// ----------------------------------------------------------------------------

// Forward for generic case.

template <typename TStream, typename TSource>
inline int
streamPut(TStream & stream, TSource const & source);

// Important special case of char.

template <typename TStream>
inline int
streamPut(Stream<TStream> & stream, char const c)
{
    return streamWriteChar(stream, c);
}

// Important special case of CharString.

template <typename TStream, typename TSpec>
inline int
streamPut(Stream<TStream> & stream, String<char, TSpec> const & source)
{
    return (streamWriteBlock(stream, toCString(source), length(source)) == length(source))  ?   0 : 1;
}

// Generic version, based on stringstream.

template <typename TStream, typename TSource>
inline int
_streamPut(Stream<TStream> & stream, TSource const & source, False const & /*tag*/)
{
    char buffer[1024] = "";
    ::std::stringstream s;

    s << source;
    if (s.fail())
        return s.fail();

    s >> buffer;
    if (s.fail())
        return s.fail();

    buffer[1023] = 0;

//TODO(h4nn3s): we should be able to use the following and then s.str() directly
// so we wouldnt need an extra buffer at all. but it doesnt work
//     s << source << std::ends;
//     if (s.fail())
//         return s.fail();

    return (streamWriteBlock(stream, buffer, strlen(buffer)) == strlen(buffer)) ? 0 : 1;
}

template <typename TStream, typename TSource>
inline int
_streamPut(TStream & target, TSource const & source, True const & /*tag*/)
{
	typename Iterator<TSource const, Standard>::Type it = begin(source, Standard());
	typename Iterator<TSource const, Standard>::Type itEnd = end(source, Standard());
    int res = 0;

	for (; it != itEnd && res == 0; ++it)
	{
		typename GetValue<TSource const>::Type val_ = getValue(it);
		res = streamPut(target, val_);
	}
	
	return res;
}

// Case: Character arrays.

template <typename TStream>
inline int
_streamPut(Stream<TStream> & stream, char const * source, True const & /*tag*/)
{
    return (streamWriteBlock(stream, source, strlen(source)) == strlen(source)) ? 0 : 1;
}

// Case: Array.
// TODO(holtgrew): Requires atEnd(it) <==> *it == 0. Remove?

template <typename TStream, typename TSourceValue>
inline int
_streamPut(TStream & stream, TSourceValue const * source, True const & /*tag*/)
{
    int res = 0;
	for (; !atEnd(source) && res == 0; ++source)
		res = _streamWrite(stream, *source);
	return res;
}

// Function entry for generic version.

template <typename TStream, typename TSource>
inline int
streamPut(TStream & stream, TSource const & source)
{
	return _streamPut(stream, source, typename IsSequence<TSource const>::Type());
}

}  // namespace seqean

#endif  // #ifndef SEQAN_STREAM_STREAM_BASE_H_
