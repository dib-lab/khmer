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
// The Single-Pass Record Reader specialization.  It works for all stream
// types.
// ==========================================================================

#ifndef SEQAN_STREAM_RECORD_READER_SINGLE_H_
#define SEQAN_STREAM_RECORD_READER_SINGLE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TFile>
class RecordReader<TFile, SinglePass<void> >;

template <typename TFile>
inline bool
_refillBuffer(RecordReader<TFile, SinglePass<void> > & recordReader);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Spec.Single-Pass RecordReader
..cat:Input/Output
..general:Class.RecordReader
..summary:Record reader specialization for single-pass reading.
..signature:RecordReader<TStream, SinglePass<void> >
..param.TStream:The @Concept.StreamConcept@ type to work on.
..remarks:Is not default or copy constructable.
..include:seqan/stream.h
 */

template <typename TFile>
class RecordReader<TFile, SinglePass<void> >
{
public:
    TFile & _file;
    unsigned _bufferSize;
    CharString _buffer;
    typedef typename Iterator<CharString, Standard>::Type TIter;
    TIter _current, _end;
    int _resultCode;
    bool _stayInOneBuffer;
    // We have to store the position of the end of the currently loaded block in the file separately because std streams
    // cannot tell their position if at end and clearing the eofbit leads to problems with atEnd().
    typedef typename Position<TFile>::Type TPosition;
    TPosition _position;  // Position in file.

    enum {
        OK = 0,
        INVALID_FORMAT
    };

    RecordReader(TFile & file)
            : _file(file), _bufferSize(BUFSIZ), _current(0), _end(0),
              _resultCode(0), _stayInOneBuffer(false), _position(0)
    {
        resize(_buffer, _bufferSize);
        _refillBuffer(*this);
    }

    RecordReader(TFile & file, unsigned bufferSize)
            : _file(file), _bufferSize(bufferSize), _current(0), _end(0),
              _resultCode(0), _stayInOneBuffer(false), _position(0)
    {
        resize(_buffer, _bufferSize);
        _refillBuffer(*this);
    }

private:
    // No default or copy constructor.
    RecordReader(RecordReader const &other): _file(other._file) {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Helper Function _refillBuffer()
// ----------------------------------------------------------------------------

template <typename TFile>
inline bool
_refillBuffer(RecordReader<TFile, SinglePass<void> > & recordReader)
{
    if (recordReader._stayInOneBuffer && recordReader._end != 0)
        // e.g. file format detection; if end==0 there hasnt yet been a buffer
        return false;
    // std::cerr << "REFILLING BUFFER" << std::endl;
    if (streamEof(recordReader._file))
        return false;
    recordReader._current = begin(recordReader._buffer, Standard());
    recordReader._position = streamTell(recordReader._file);
    size_t bytesRead = streamReadBlock(recordReader._current, recordReader._file, recordReader._bufferSize);
    recordReader._position += bytesRead;
    if (bytesRead != recordReader._bufferSize) {
        // If we read fewer characters and the stream is not at its end then
        // there was an error reading the file.
        recordReader._resultCode = streamError(recordReader._file);
        if (recordReader._resultCode) {
            // std::cerr << "RESULT IS " << recordReader._resultCode << std::endl;
            recordReader._end = recordReader._current;
            return false;
        }
    }
    recordReader._end = recordReader._current + bytesRead;
    // std::cerr << "read " << bytesRead << " bytes" << std::endl;
    // std::cerr << "ERROR? " << streamError(recordReader._file) << std::endl;
    // std::cerr << "EOF? " << streamEof(recordReader._file) << std::endl;

    return true;
}

// ----------------------------------------------------------------------------
// Function position()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Document!
template <typename TFile>
inline typename Position<TFile>::Type
position(RecordReader<TFile, SinglePass<void> > const & recordReader)
{
    typename Position<TFile>::Type bufferedUnread = recordReader._end - recordReader._current;
    return recordReader._position - bufferedUnread;
}

// ----------------------------------------------------------------------------
// Function setPosition()
// ----------------------------------------------------------------------------

template <typename TFile, typename TPosition>
inline int
setPosition(RecordReader<TFile, SinglePass<void> > & recordReader, TPosition pos)
{
    int res = streamSeek(recordReader._file, pos, SEEK_SET);
    if (res != 0)
        return res;
    _refillBuffer(recordReader);
    return 0;
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

template <typename TFile>
inline bool
atEnd(RecordReader<TFile, SinglePass<void> > & recordReader)
{
    // There is more data if the buffer is not exhausted.
    if (recordReader._current != recordReader._end)
        return false;
    // There is no more data if the buffer is exhausted and there the stream
    // is at the end of the file or there previously was an error reading the
    // file.
    
    if (streamEof(recordReader._file) || recordReader._resultCode != 0)
        return true;

    // std::cerr << "refilling in atEnd()" << std::endl;

    // Otherwise, we can try to load some data.  This case only happens if
    // atEnd is called at the beginning of the file; Otherwise, goNext()
    // will load data.
    return !_refillBuffer(recordReader);
}

// ----------------------------------------------------------------------------
// Function resultCode()
// ----------------------------------------------------------------------------

template <typename TFile>
inline int
resultCode(RecordReader<TFile, SinglePass<void> > & recordReader)
{
    return recordReader._resultCode;
}

// ----------------------------------------------------------------------------
// Function goNext()
// ----------------------------------------------------------------------------

template <typename TFile>
inline bool
goNext(RecordReader<TFile, SinglePass<void> > & recordReader)
{
    SEQAN_ASSERT(recordReader._current != recordReader._end);

    recordReader._current += 1;
    // If there is more data in the buffer then we're done.
    if (recordReader._current != recordReader._end)
        return false;  // Has more data.

    // Otherwise, try to load some data.
    return !_refillBuffer(recordReader);
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TFile>
inline char
value(RecordReader<TFile, SinglePass<void> > & recordReader)
{
    SEQAN_ASSERT(recordReader._current != recordReader._end);
    return *recordReader._current;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_RECORD_READER_SINGLE_H_
