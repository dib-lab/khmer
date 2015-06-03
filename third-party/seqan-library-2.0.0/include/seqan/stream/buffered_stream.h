// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_STREAM_BUFFERED_STREAM_
#define SEQAN_STREAM_BUFFERED_STREAM_

#include <cstdio>
#include <cstring>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TValue, typename TDirection, typename TTraits = std::char_traits<TValue> >
class BufferedStreamBuf;

// ============================================================================
// Tags, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class BufferedStream
// ----------------------------------------------------------------------------

// TODO(esiragusa): @extends BasicStream::Type

/*!
 * @class BufferedStream
 * @implements StreamConcept
 * @headerfile <seqan/stream.h>
 * @brief Adds a buffer to another @link StreamConcept stream @endlink.
 *
 * @signature template <typename TUnbufferedStream, TDirection>
 *            class BufferedStream;
 *
 * @tparam TUnbufferedStream The type of the unbuffered @link StreamConcept stream @endlink to wrap.
 * @tparam TDirection        The stream direction, one of @link DirectionTags @endlink.
 */

template <typename TUnbufferedStream, typename TDirection>
class BufferedStream :
    public BasicStream<
        typename TUnbufferedStream::char_type,
        TDirection,
        typename TUnbufferedStream::traits_type>::Type
{
protected:
    typedef typename TUnbufferedStream::char_type                   TValue;
    typedef typename TUnbufferedStream::traits_type                 TTraits;
    typedef typename BasicStream<TValue, TDirection, TTraits>::Type TBasicStream;
    typedef std::basic_streambuf<TValue, TTraits>                   TStreamBuf;

    BufferedStreamBuf<TValue, TDirection, TTraits> buf;

public:

    /*!
     * @fn BufferedStream::BufferedStream
     * @brief Default constructor and construction from to-be-wrapped unbuffered stream.
     *
     * @signature BufferedStream::BufferedStream();
     * @signature BufferedStream::BufferedStream(unbufferedStream);
     *
     * @param[in] stream The to-be-wrapped stream.
     */
    BufferedStream() :
        TBasicStream(&buf)
    {}

    explicit BufferedStream(TUnbufferedStream & stream) :
        TBasicStream(&buf), buf(stream.rdbuf())
    {}

    /*!
     * @fn BufferedStream::setStreamBuf
     * @brief Set the stream buffer of the BufferedStream.
     *
     * @signature void BufferedStream::setStreamBuf(streamBuf);
     *
     * @param[in] streamBuf The <tt>std::basic_streambuf&lt;&gt;</tt> to use.
     */
    void setStreamBuf(TStreamBuf & streamBuf)
    {
        buf.setStreamBuf(streamBuf);
    }

    /*!
     * @fn BufferedStream::setStream
     * @brief Set the underlying stream to wrap.
     *
     * @signature void BufferedStream::setStream(stream);
     *
     * @param[in] stream The <tt>TUnbufferedStream</tt> to use.
     */
    void setStream(TUnbufferedStream & stream)
    {
        setStreamBuf(*stream.rdbuf());
    }
};

// ----------------------------------------------------------------------------
// Class BufferedStream
// ----------------------------------------------------------------------------

// TODO(holtgrew): Implementation detail, should thus be called BufferedStreamBuf_ or documented.

template <typename TValue, typename TDirection, typename TTraits>
class BufferedStreamBuf :
    public std::basic_streambuf<TValue, TTraits>
{
protected:
    typedef std::basic_streambuf<TValue, TTraits>                   TStreamBuf;
    typedef typename TTraits::int_type                              TInt;

    static const size_t defaultBufferSize = 1024;   // size of the data buffer
    static const size_t defaultPutbackSize = 256;   // size of the buffer that should be used as putback area
    std::vector<TValue> buffer;                     // data buffer
    size_t bufferSize;
    size_t putbackSize;

    using TStreamBuf::eback;
    using TStreamBuf::gptr;
    using TStreamBuf::egptr;
    using TStreamBuf::setg;

public:
    TStreamBuf *streamBufPtr;

    /* constructor
     * - initialize empty data buffer
     * - no putback area
     * => force underflow()
     */
    BufferedStreamBuf(size_t bufferSize = defaultBufferSize,
                      size_t putbackSize = defaultPutbackSize) :
        bufferSize(bufferSize),
        putbackSize(putbackSize),
        streamBufPtr(NULL)
    {
        clear();
    }

    explicit
    BufferedStreamBuf(TStreamBuf &streamBuf,
                      size_t bufferSize = defaultBufferSize,
                      size_t putbackSize = defaultPutbackSize) :
        bufferSize(bufferSize),
        putbackSize(putbackSize)
    {
        setStreamBuf(streamBuf);
    }

    void clear()
    {
        setg((TValue*)NULL,     // beginning of putback area
             (TValue*)NULL,     // read position
             (TValue*)NULL);    // end position
    }

    void setStreamBuf(TStreamBuf &streamBuf)
    {
        clear();
        streamBufPtr = &streamBuf;
    }

protected:
    // insert new characters into the buffer
    virtual TInt underflow()
    {
        // is read position before end of buffer?
        if (gptr() < egptr())
            return TTraits::to_int_type(*gptr());

        /* process size of putback area
         * - use number of characters read
         * - but at most four
         */
        size_t numPutback;
        numPutback = gptr() - eback();
        if (numPutback > putbackSize)
            numPutback = putbackSize;

        buffer.resize(bufferSize);

        /* copy up to four characters previously read into
         * the putback buffer (area of first four characters)
         */
        std::memmove(
            &buffer[putbackSize - numPutback],
            gptr() - numPutback,
            numPutback * sizeof(TValue));

        // read new characters
        size_t numRead = streamBufPtr->sgetn(
                            &buffer[putbackSize],
                            buffer.size() - putbackSize);

        // reset buffer pointers
        setg(&buffer[putbackSize - numPutback],     // beginning of putback area
             &buffer[putbackSize],                  // read position
             &buffer[putbackSize + numRead]);       // end position

        if (numRead <= 0)
        {
            // ERROR or EOF
            return EOF;
        }

        // return next character
        return TTraits::to_int_type(*gptr());
    }
};

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_BUFFERED_STREAM_
