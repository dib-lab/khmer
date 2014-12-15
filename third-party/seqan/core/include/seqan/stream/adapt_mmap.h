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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// Adaptions for the Memorymapped Strings
// ==========================================================================

// TODO(holtgrew): Should better be string adaption!

#ifndef SEQAN_STREAM_ADAPT_MMAP_H_
#define SEQAN_STREAM_ADAPT_MMAP_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*
.Adaption.IO stream
..cat:Input/Output
..summary:Adaption from $fstream$, $ifstream$ and $ofstream$ to the @Concept.StreamConcept@ concept.
..include:seqan/stream.h
 */

// ============================================================================
// Metafunctions
// ============================================================================


// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function streamWriteChar()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TChar>
inline int
streamWriteChar(String<TValue, TSpec> & stream, TChar const & c)
{
    appendValue(stream, c);
    return 0;
}

// ----------------------------------------------------------------------------
// Function streamWriteBlock()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline typename Size<String<TValue, TSpec> >::Type
streamWriteBlock(String<TValue, TSpec> & stream, char const * ptr, unsigned count)
{
    reserve(stream, length(stream) + count);
    for (unsigned i = 0; i < count; ++i, ++ptr)
        appendValue(stream, *ptr);
    return count;
}

// ----------------------------------------------------------------------------
// Function streamPut()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline int
streamPut(String<TValue, TSpec> & stream, char const c)
{
    appendValue(stream, c);
    return 0;
}

template <typename TValue1, typename TSpec1, typename TValue2, typename TSpec2>
inline int
streamPut(String<TValue1, TSpec1> & stream,
          SimpleType<TValue2, TSpec2> const & c)
{
    appendValue(stream, c);
    return 0;
}


// template <typename TValue0, typename TSpec0,
//           typename TValue, typename TSpec, typename TSpec2>
// inline int
// streamPut(String<TValue0, MMap<TSpec0> > & stream,
//           String<SimpleType<TValue, TSpec>, TSpec2> const & source)
// {
//     String<char, CStyle> buf = source;
//     append(stream, toCString(buf));
//     return 0;
// }

template <typename TValue, typename TSpec, typename TValue2, typename TSpec2>
inline int
streamPut(String<TValue, TSpec> & stream,
          String<TValue2, TSpec2> const & source)
{
    append(stream, source);
    return 0;
}


template <typename TValue, typename TSpec, typename TSource>
inline int
_appendWithoutTrailing0(String<TValue, TSpec> & stream,
                        TSource const & source)
{
    for (int i = 0; source[i] != 0; ++i)
        appendValue(stream, source[i]);
    return 0;
}

template <typename TValue, typename TSpec>
inline int
streamPut(String<TValue, TSpec> & stream, char const *source)
{
    return _appendWithoutTrailing0(stream, source);
}


// for numerical types
template <typename TValue, typename TSpec, typename TSource>
inline int
streamPut(String<TValue, TSpec> & stream, TSource const & source)
{
    std::ostringstream str;
    str << source << std::ends;
    if (str.fail())
        return str.fail();
    return _appendWithoutTrailing0(stream, str.str());
}

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_ADAPT_MMAP_H_
