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
// Interface for the stream concept.
// ==========================================================================

#ifndef SEQAN_STREAM_CONCEPT_STREAM_H_
#define SEQAN_STREAM_CONCEPT_STREAM_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Concept.StreamConcept
..summary:Concept for I/O streams.
 */

// ============================================================================
// Metafunctions
// ============================================================================

/**
.Tag.Stream Feature
..cat:Input/Output
..summary:Tag to select a given feature for querying in @Metafunction.HasStreamFeature@.
..tag.IsInput:Query whether the stream supports @Function.streamReadChar@ and @Function.streamReadBlock@.
..tag.IsOutput:Query whether the stream supports @Function.streamWriteChar@ and @Function.streamWriteBlock@.
..tag.HasPeek:Query whether the stream supports @Function.streamPeek@.
..tag.HasFilename:Query whether the stream has a file name.
..tag.Seek<OriginBegin>:Query whether it is possible to @Function.streamSeek@ with $SEEK_SET$.
..tag.Seek<OriginCurrent>:Query whether it is possible to @Function.streamSeek@ with $SEEK_CUR$.
..tag.Seek<OriginEnd>:Query whether it is possible to @Function.streamSeek@ with $SEEK_END.
..tag.Tell:Query whether the stream has a @Function.streamTell@ function.
..include:seqan/stream.h
 */

struct IsInput_;
typedef Tag<IsInput_> IsInput;

struct IsOutput_;
typedef Tag<IsOutput_> IsOutput;

struct HasPeek_;
typedef Tag<HasPeek_> HasPeek;

struct HasFilename_;
typedef Tag<HasFilename_> HasFilename;

struct OriginBegin_;
typedef Tag<OriginBegin_> OriginBegin;

struct OriginEnd_;
typedef Tag<OriginEnd_> OriginEnd;

struct OriginCurrent_;
typedef Tag<OriginCurrent_> OriginCurrent;

template <typename TSpec>
struct Seek;

struct Tell_;
typedef Tag<Tell_> Tell;

/**
.Metafunction.HasStreamFeature
..concept:Concept.StreamConcept
..cat:Input/Output
..summary:Query features of a stream type.
..signature:HasStreamFeature<TStream, TFeatureTag>::Type
..param.TStream:The stream type to query the property of.
...type:Concept.StreamConcept
..param.TFeatureTag:The feature tag.
...type:Tag.Stream Feature
..returns:Either $True$ or $False$.
..remarks:Note that this only checks whether the type principally has this feature. For example, if a stream wraps a Unix PIPE internally, it might not support seek.
..include:seqan/stream.h
 */

template <typename TStream, typename TTag>
struct HasStreamFeature;

// ============================================================================
// Functions
// ============================================================================

/**
.Function.streamPeek
..concept:Concept.StreamConcept
..cat:Input/Output
..summary:Read next character from stream without advancing current position.
..signature:streamPeek(c, stream)
..param.c:The read character is written here.
...type:nolink:$char &$
..param.stream:The stream object to read from.
...type:Concept.StreamConcept
..returns:$int$, 0 on success, otherwise the error value from the underlying string system.
..remarks:Note that this might involve two calls into the stream library, e.g. for cstdio streams, it involves a call to both $getc()$ and $ungetc()$.
..see:Function.streamReadChar
..see:Function.streamReadBlock
..include:seqan/stream.h

.Function.streamReadChar
..concept:Concept.StreamConcept
..cat:Input/Output
..summary:Read next character from stream and advance the current position.
..signature:streamReadChar(c, stream)
..param.c:The read character is written here.
...type:nolink:$char &$
..param.stream:The stream object to read from.
...type:Concept.StreamConcept
..returns:$int$, 0 on success, otherwise the error value from the underlying string system.
..see:Function.streamPeek
..see:Function.streamReadBlock
..include:seqan/stream.h

.Function.streamEof
..concept:Concept.StreamConcept
..cat:Input/Output
..summary:Check end-of-file state of a @Concept.StreamConcept@.
..signature:streamReadChar(stream)
..param.stream:The stream object to read from.
...type:Concept.StreamConcept
..returns:$bool$, true if the stream is in end-of-file state, false otherwise.
..remarks:Note that the exact behaviour depends on the underlying implementation. With $FILE *$ streams, the stream can be in the end-of-file state from the point where the last character was read or from the point where the user tried to read beyond the end of the file.
..include:seqan/stream.h

.Function.streamReadBlock
..concept:Concept.StreamConcept
..cat:Input/Output
..summary:Read a block of bytes into a buffer.
..signature:streamReadBlock(target, stream, maxLen)
..param.target:The buffer to read into. There has to be enough space for $maxLen$ bytes.
...type:nolink:$char *$
..param.stream:The stream to read from.
...type:Concept.StreamConcept
..param.maxLen:maximal number of characters to read.
...type:nolink:$size_t$
..returns:Number of read bytes.
..example.text:Copying data from a std::fstream into another std::fstream using SeqAn's stream adaption.
..example.code:
#include <fstream>
#include <seqan/sequence.h>
#include <seqan/stream.h>

int main()
{
    std::fstream in("in.txt", std::ios::binary | std::ios::in);
    std::fstream out("out.txt", std::ios::binary | std::ios::in);

    seqan::CharString buffer;
    resize(buffer, 1000);

    while (!seqan::atEnd(in) && seqan::streamError(in) == 0)
    {
        int num = seqan::streamReadBlock(&buffer[0], in, length(buffer));
        seqan::streamWriteBlock(out, &buffer[0], num);
    }

    return 0;
}
..see:Function.streamPeek
..see:Function.streamReadChar
..include:seqan/stream.h

.Function.streamWriteChar
..concept:Concept.StreamConcept
..cat:Input/Output
..summary:Write one character to the stream.
..signature:streamWriteChar(stream, c)
..param.stream:The stream object to write to.
...type:Concept.StreamConcept
..param.c:The character to write to the stream.
...type:nolink:$char$
..returns:$int$ with error code, 0 on success.
..see:Function.streamWriteBlock
..include:seqan/stream.h

.Function.streamWriteBlock
..concept:Concept.StreamConcept
..cat:Input/Output
..summary:Write a block of bytes from a buffer into a stream.
..signature:streamWriteBlock(stream, source, count)
..param.stream:The stream object to write to.
...type:Concept.StreamConcept
..param.source:The data to write to the stream.
...type:nolink:$char *$
..param.count:The number of bytes to write to the stream.
...type:nolink:$size_t$
..returns:$int$ with the number of successfully written objects. 
..example.text:Copying data from a std::fstream into another std::fstream using SeqAn's stream adaption.
..example.code:
#include <fstream>
#include <seqan/sequence.h>
#include <seqan/stream.h>

int main()
{
    std::fstream in("in.txt", std::ios::binary | std::ios::in);
    std::fstream out("out.txt", std::ios::binary | std::ios::in);

    seqan::CharString buffer;
    resize(buffer, 1000);

    while (!seqan::atEnd(in) && seqan::streamError(in) == 0)
    {
        int num = seqan::streamReadBlock(&buffer[0], in, length(buffer));
        seqan::streamWriteBlock(out, &buffer[0], num);
    }

    return 0;
}
..see:Function.streamWriteChar
..include:seqan/stream.h

.Function.streamPut
..concept:Concept.StreamConcept
..cat:Input/Output
..summary:Write different types to stream
..signature:streamPut(stream, source)
..param.stream:The stream object to write to.
...type:Concept.StreamConcept
..param.source:The data to write to the stream.
...type:nolink:$char *$ (must be 0-terminated!)
...type:Shortcut.CharString
...type:nolink:String<char, *>
...type:nolink:numerical types ($int$, $double$ ...)
..returns:$int$ with an error code, 0 on success.
..remarks:Implementation note: for some specializations of @Concept.StreamConcept@ certain conversions take place through stringstream and a buffer of size 1023. It follows that the result of the conversion cannot be longer. However this should only effect numericals right now. If you still encounter truncated strings with another type, convert to $const char*$ manually before writing.
..see:Function.streamWriteChar
..see:Function.streamWriteBlock
..include:seqan/stream.h

.Function.streamFlush
..concept:Concept.StreamConcept
..cat:Input/Output
..summary:Flush the underlying stream.
..signature:streamFlush(stream)
..param.stream:The stream object to flush.
...type:Concept.StreamConcept
..returns:$int$ with an error code, 0 on success.
..include:seqan/stream.h

.Function.streamError
..concept:Concept.StreamConcept
..cat:Input/Output
..summary:Return the stream's error code.
..signature:streamError(stream)
..param.stream:The stream object to query.
...type:Concept.StreamConcept
..returns:$int$ with an error code, 0 is used for "no errors."
..include:seqan/stream.h

.Function.streamSeek
..concept:Concept.StreamConcept
..cat:Input/Output
..summary:Perform a seek operation on the stream.
..signature:streamSeek(stream, delta, origin)
..param.stream:The stream object to seek on.
...type:Concept.StreamConcept
..param.delta:The (relative) position.
...type:nolink:$long int$
..param.origin:Where to query from, we use the integers from <cstdio>.
...type:nolink:$int$
...remarks:Possible values are $SEEK_SET$, $SEEK_CUR$, and $SEEK_END$.
..returns:$int$ with an error code, 0 is used for "no errors."
..see:Function.streamTell
..include:seqan/stream.h

.Function.streamTell
..concept:Concept.StreamConcept
..cat:Input/Output
..summary:Get the position in the current stream.
..signature:streamTell(stream)
..param.stream:The stream object to query the current position for.
...type:Concept.StreamConcept
..returns:The position within the stream, of type @Metafunction.Position@.
..see:Function.streamSeek
..include:seqan/stream.h
 */

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_CONCEPT_STREAM_H_
