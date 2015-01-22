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
// The Single-Pass Record Reader working on strings.  It uses that data does
// not have to be read explicitely but is automatically loaded by the OS on
// pagefaults.
// ==========================================================================

#ifndef SEQAN_STREAM_RECORD_READER_SINGLE_MMAP_H_
#define SEQAN_STREAM_RECORD_READER_SINGLE_MMAP_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// TODO(holtgrew): Also use mmadvise to mark used memory as such?

// TODO(holtgrew): This could easily be adjusted to work for any string specialization by adding another layer, signature would then be RecordReader<TFile, SinglePass<StringReader<TStringSpec> > >.

/**
.Spec.Single-Pass String RecordReader
..cat:Input/Output
..general:Spec.Single-Pass RecordReader
..summary:Record reader specialization for single-pass reading from memory mapped files.
..signature:RecordReader<TStream, SinglePass<StringReader> >
..param.TStream:The @Concept.StreamConcept@ type to work on.
..remarks:This record reader does not have any buffers but uses the memory mapped string directly.
..remarks:Is not default or copy constructable.
..remarks:The buffer size is the granularity in which @Function.mmapAdvise@ will be called if $TString$ is an @Spec.MMap String@.
..see:Class.String
..include:seqan/stream.h
 */

// TODO(holtgrew): Actually, this should work on sequences but is yet untested.

template <typename TString>
class RecordReader<TString, SinglePass<StringReader> >
{
public:
    typedef typename Iterator<TString, Standard>::Type TIter;
    typedef typename Size<TString>::Type TSize;

    TString & _string;
    TIter _current, _end;
    TSize _bufferSize;

    enum {
        OK = 0,
        INVALID_FORMAT
    };

    RecordReader(TString & string)
            : _string(string), _current(begin(string)), _end(end(string)),
              _bufferSize(BUFSIZ)
    {}

    RecordReader(TString & string, unsigned bufferSize)
            : _string(string), _current(begin(string)), _end(end(string)),
              _bufferSize(bufferSize)
    {}

private:
    // No default or copy constructor.
    RecordReader(RecordReader const &other): _string(other._string) {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function position()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Document!
template <typename TString>
inline typename Position<TString>::Type
position(RecordReader<TString, SinglePass<StringReader> > const & recordReader)
{
    return recordReader._current - begin(recordReader._string, Standard());
}

// ----------------------------------------------------------------------------
// Function setPosition()
// ----------------------------------------------------------------------------

template <typename TString, typename TPosition>
inline int
setPosition(RecordReader<TString, SinglePass<StringReader> > & recordReader, TPosition pos)
{
    recordReader._current = iter(recordReader._string, pos, Standard());
    return 0;
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

template <typename TString>
inline bool
atEnd(RecordReader<TString, SinglePass<StringReader> > & recordReader)
{
    return recordReader._current == recordReader._end;
}

// TODO(holtgrew): We would rather have the const version only.

template <typename TString>
inline bool
atEnd(RecordReader<TString, SinglePass<StringReader> > const & recordReader)
{
    return recordReader._current == recordReader._end;
}

// ----------------------------------------------------------------------------
// Function resultCode()
// ----------------------------------------------------------------------------

template <typename TString>
inline int
resultCode(RecordReader<TString, SinglePass<StringReader> > & /*recordReader*/)
{
    return 0;
}

// ----------------------------------------------------------------------------
// Function goNext()
// ----------------------------------------------------------------------------

template <typename TString>
inline bool
goNext(RecordReader<TString, SinglePass<StringReader> > & recordReader)
{
    SEQAN_ASSERT(recordReader._current != recordReader._end);
    recordReader._current += 1;
    return false;
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TString>
inline char
value(RecordReader<TString, SinglePass<StringReader> > & recordReader)
{
    SEQAN_ASSERT(recordReader._current != recordReader._end);
    return *recordReader._current;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_RECORD_READER_SINGLE_MMAP_H_
