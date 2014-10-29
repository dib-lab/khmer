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
// The Double-Pass Record Reader specialization.  It works for all stream
// types.
// ==========================================================================

#ifndef SEQAN_STREAM_RECORD_READER_DOUBLE_H_
#define SEQAN_STREAM_RECORD_READER_DOUBLE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Spec.Double-Pass RecordReader
..cat:Input/Output
..general:Class.RecordReader
..summary:Record reader specialization for double-pass reading.
..signature:RecordReader<TStream, DoublePass>
..param.TStream:The @Concept.StreamConcept@ type to work on.
..remarks:Is not default or copy constructable.
..include:seqan/stream.h
 */

template <typename TFile>
class RecordReader<TFile, DoublePass<> >
{
public:
    TFile & _file;
    unsigned _bufferSize;
    // INVARIANT: No used buffer may be empty after reading into it!
    String<CharString *> _usedBuffers;    // Buffers with relevant data.
    // TODO(holtgrew): It is probably a bad idea to do our own pooling here!
    String<CharString *> _unusedBuffers;  // Buffers without relevant data.
    typedef typename Iterator<CharString, Standard>::Type TIter;
    TIter _current, _end;
    CharString * _currentBuffer;
    unsigned _currentBuffNo;
    int _resultCode;
    int _passNo;
    char * _beginInFirst;
    bool _stayInOneBuffer; // needed for stream format detection
    // In contrast to the SinglePass<> reader, we store the current position of the file in _position since it is more
    // complex to do the computation otherwise.  This means more updates to _position, though.
    typedef typename Position<TFile>::Type TPosition;
    TPosition _position;  // Position in file.
    TPosition _firstPassPos;  // Position of _firstPass call.

    enum {
        OK = 0,
        INVALID_FORMAT
    };

    RecordReader(TFile & file)
            : _file(file), _bufferSize(BUFSIZ), _current(0), _end(0),
              _currentBuffer(0), _currentBuffNo(0), _resultCode(0),
              _passNo(0), _beginInFirst(0), _stayInOneBuffer(false),
              _position(0), _firstPassPos(0)
    {
        // resize(_buffer, _bufferSize);
    }

    RecordReader(TFile & file, unsigned bufferSize)
            : _file(file), _bufferSize(bufferSize), _current(0), _end(0),
              _currentBuffer(0), _currentBuffNo(0), _resultCode(0),
              _passNo(0), _beginInFirst(0), _stayInOneBuffer(false),
              _position(0), _firstPassPos(0)
    {
        // resize(_buffer, _bufferSize);
    }

    ~RecordReader()
    {
        for (unsigned i = 0; i < length(_usedBuffers); ++i)
            delete _usedBuffers[i];
        for (unsigned i = 0; i < length(_unusedBuffers); ++i)
            delete _unusedBuffers[i];
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
// Helper Function _fillNextBuffer()
// ----------------------------------------------------------------------------

template <typename TFile>
inline bool
_fillNextBuffer(RecordReader<TFile, DoublePass<> > & recordReader)
{
    if (recordReader._stayInOneBuffer && recordReader._end != 0)
    // e.g. FileFormat Detection
        return false;
    // std::cerr << "REFILLING BUFFER" << std::endl;
    if (streamEof(recordReader._file))
        return false;
    if (empty(recordReader._unusedBuffers)) {
        recordReader._currentBuffer = new CharString();
        appendValue(recordReader._usedBuffers, recordReader._currentBuffer);
    } else {
        recordReader._currentBuffer = back(recordReader._unusedBuffers);
        eraseBack(recordReader._unusedBuffers);
        appendValue(recordReader._usedBuffers, recordReader._currentBuffer);
    }

    resize(*recordReader._currentBuffer, recordReader._bufferSize);
    recordReader._current = begin(*recordReader._currentBuffer, Standard());
    // std::cerr << "recordReader._current = begin(...) == " << (void*)(recordReader._current) << " [fill next buffer]" << std::endl;
    recordReader._currentBuffNo += 1;

    if (recordReader._beginInFirst == NULL)
        recordReader._beginInFirst = recordReader._current;

    size_t bytesRead = streamReadBlock(recordReader._current, recordReader._file, recordReader._bufferSize);
    if (bytesRead != recordReader._bufferSize) {
        resize(*recordReader._currentBuffer, bytesRead);
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
// Function startFirstPass()
// ----------------------------------------------------------------------------

/**
.Function.startFirstPass
..class:Spec.Double-Pass RecordReader
..cat:Input/Output
..signature:startFirstPass(recordReader)
..summary:Start the first reading pass.
..param.recordReader:RecordReader to start the first pass with.
...type:Spec.Double-Pass RecordReader
..remarks:This will memoize the current position in the buffer, to return to in second pass.
..see:Function.startSecondPass
..include:seqan/stream.h
 */

template <typename TFile>
void
startFirstPass(RecordReader<TFile, DoublePass<> > & recordReader)
{
    // Store begin position of first pass.
    recordReader._firstPassPos = recordReader._position;
    recordReader._passNo = 1;
    recordReader._beginInFirst = recordReader._current;
    // TODO(holtgrew_): Add assertion that previous second pass ended at same position.
    // What if is exactly at end of buffer?
    if (length(recordReader._usedBuffers) > 1u) {
        CharString * ptr = back(recordReader._usedBuffers);
        eraseBack(recordReader._usedBuffers);
        append(recordReader._unusedBuffers, recordReader._usedBuffers);
        clear(recordReader._usedBuffers);
        appendValue(recordReader._usedBuffers, ptr);
    }

    // Fill the next (or even: first) buffer if the current is empty.
    if (recordReader._current == recordReader._end)
        _fillNextBuffer(recordReader);
    // std::cerr << "starting first pass from " << (void*)(recordReader._current) << std::endl;
}

// ----------------------------------------------------------------------------
// Helper Function _jumpToNextBuffer()
// ----------------------------------------------------------------------------

template <typename TFile>
bool
_jumpToNextBuffer(RecordReader<TFile, DoublePass<> > & recordReader)
{
    if (recordReader._stayInOneBuffer) // e.g. file format detection
        return false;
    SEQAN_ASSERT_EQ(recordReader._passNo, 2);
    SEQAN_ASSERT(recordReader._current ==  recordReader._end);
    if (recordReader._currentBuffNo + 1 >= length(recordReader._usedBuffers))
        return false;
    recordReader._currentBuffNo += 1;
    recordReader._currentBuffer = recordReader._usedBuffers[recordReader._currentBuffNo];
    recordReader._current = begin(*recordReader._currentBuffer, Standard());
    // std::cerr << "recordReader._current = begin(...) == " << (void*)(recordReader._current) << " [jump to next]" << std::endl;
    recordReader._end = end(*recordReader._currentBuffer, Standard());
    return true;
}

// ----------------------------------------------------------------------------
// Function startSecondPass()
// ----------------------------------------------------------------------------

/**
.Function.startSecondPass
..class:Spec.Double-Pass RecordReader
..cat:Input/Output
..signature:startSecondPass(recordReader)
..summary:Start the second reading pass.
..param.recordReader:RecordReader to start the second pass with.
...type:Spec.Double-Pass RecordReader
..remarks:This will reset the position in the buffer
..see:Function.startFirstPass
..include:seqan/stream.h
 */

template <typename TFile>
void
startSecondPass(RecordReader<TFile, DoublePass<> > & recordReader)
{
    // Set position back to first pass start position.
    recordReader._position = recordReader._firstPassPos;
    
    SEQAN_ASSERT_EQ(recordReader._passNo, 1);
    recordReader._passNo = 2;
    recordReader._currentBuffNo = 0;
    recordReader._currentBuffer = recordReader._usedBuffers[0];
    recordReader._current = recordReader._beginInFirst;
    // std::cerr << "recordReader._current = " << (void*)(recordReader._beginInFirst) << " [start second pass]" << std::endl;
    recordReader._end = end(*recordReader._currentBuffer, Standard());
}

// ----------------------------------------------------------------------------
// Function position()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Document!
template <typename TFile>
inline typename Position<TFile>::Type
position(RecordReader<TFile, DoublePass<void> > const & recordReader)
{
    return recordReader._position;
}

// ----------------------------------------------------------------------------
// Function setPosition()
// ----------------------------------------------------------------------------

// This automatically starts the first pass at the position.
template <typename TFile, typename TPosition>
inline int
setPosition(RecordReader<TFile, DoublePass<void> > & recordReader, TPosition pos)
{
    // Clear buffers, mark as unused.
    append(recordReader._unusedBuffers, recordReader._usedBuffers);
    clear(recordReader._usedBuffers);

    // Seek to position in file.
    int res = streamSeek(recordReader._file, pos, SEEK_SET);
    if (res != 0)
        return res;
    recordReader._position = pos;
    recordReader._current = recordReader._end = 0;
    startFirstPass(recordReader);
    return 0;
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

template <typename TFile>
inline bool
atEnd(RecordReader<TFile, DoublePass<> > & recordReader)
{
    // std::cerr << "HAS MORE " << __LINE__ << std::endl;
    // There can be more data if the current buffer is not exhausted.
    if (recordReader._current != recordReader._end)
        return false;
    // std::cerr << "HAS MORE " << __LINE__ << std::endl;
    // There can be more data if the current buffer is not the last one.
    if (recordReader._currentBuffNo + 1 < length(recordReader._usedBuffers))
        return false;
    // std::cerr << "HAS MORE " << __LINE__ << std::endl;

    // There is no more data if the buffer is exhausted and the stream
    // is at the end of the file or there previously was an error reading the
    // file.

    if (streamEof(recordReader._file) || recordReader._resultCode != 0)
        return true;
    // std::cerr << "HAS MORE " << __LINE__ << std::endl;

    // std::cerr << "refilling in atEnd()" << std::endl;

    // Otherwise, we can try to load some data.  This case only happens if
    // atEnd() is called at the beginning of the file; Otherwise, goNext()
    // will load data.
    switch (recordReader._passNo) {
        case 0:
            return false;
        case 1:
            return !_fillNextBuffer(recordReader);
        case 2:
            if (!_jumpToNextBuffer(recordReader))
                return !_fillNextBuffer(recordReader);
            else
                return false;
        default:
            SEQAN_ASSERT_FAIL("Invalid pass no: %d", recordReader._passNo);
            return true;
    }
}

// ----------------------------------------------------------------------------
// Function resultCode()
// ----------------------------------------------------------------------------

template <typename TFile>
inline int
resultCode(RecordReader<TFile, DoublePass<> > & recordReader)
{
    return recordReader._resultCode;
}

// ----------------------------------------------------------------------------
// Function goNext()
// ----------------------------------------------------------------------------

template <typename TFile>
inline bool
goNext(RecordReader<TFile, DoublePass<> > & recordReader)
{
    SEQAN_ASSERT(recordReader._current != recordReader._end);

    // std::cerr << "recordReader._current += 1 == " << (void*)(recordReader._current) << std::endl;
    recordReader._current += 1;
    recordReader._position += 1;
    // If there is more data in the buffer then we're done.
    if (recordReader._current != recordReader._end)
        return false;  // Has more data.

    // std::cerr << "REFILLING IN goNext()" << std::endl;

    SEQAN_ASSERT_GEQ(recordReader._passNo, 1);
    SEQAN_ASSERT_LEQ(recordReader._passNo, 2);
    if (recordReader._passNo == 1)
        return !_fillNextBuffer(recordReader);
    else
        return !_jumpToNextBuffer(recordReader);
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TFile>
inline char
value(RecordReader<TFile, DoublePass<> > & recordReader)
{
    SEQAN_ASSERT(recordReader._current != recordReader._end);
    // std::cerr << "*recordReader._current == " << *recordReader._current << std::endl;
    return *recordReader._current;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_RECORD_READER_DOUBLE_H_
