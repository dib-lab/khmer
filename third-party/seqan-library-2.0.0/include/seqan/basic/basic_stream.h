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
// Basic definitions for the stream module.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_BASIC_STREAM_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_BASIC_STREAM_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TDirection>
struct StreamIterator;

template <typename TValue, typename TTraits, typename TValue2>
inline void writeValue(std::basic_ostream<TValue, TTraits> &ostream, TValue2 val);

template <typename TValue, typename TTraits, typename TValue2>
inline void writeValue(std::ostreambuf_iterator<TValue, TTraits> &iter, TValue2 val);

template <typename TContainer, typename TValue>
inline void writeValue(Iter<TContainer, StreamIterator<Output> > &iter, TValue val);

template <typename TValue, typename TTraits>
inline bool atEnd(std::istreambuf_iterator<TValue, TTraits> const &it);

template <typename TChar, typename TCharTraits, typename TAlloc>
inline typename std::basic_string<TChar, TCharTraits, TAlloc>::size_type
length(std::basic_string<TChar, TCharTraits, TAlloc> const & me);

/*!
 * @macro SEQAN_HAS_ZLIB
 * @headerfile <seqan/stream.h>
 * @brief Defined as 0 or 1, depending on zlib being available.
 *
 * @signature #define SEQAN_HAS_ZLIB 0  // or 1
 */

/*!
 * @macro SEQAN_HAS_BZIP2
 * @headerfile <seqan/stream.h>
 * @brief Defined as 0 or 1, depending on bzlib being available.
 *
 * @signature #define SEQAN_HAS_BZIP 0  // or 1
 */

// ============================================================================
// Tags
// ============================================================================

// ============================================================================
// Concepts
// ============================================================================

// --------------------------------------------------------------------------
// Concept StreamConcept
// --------------------------------------------------------------------------

/*!
 * @concept StreamConcept
 * @headerfile <seqan/basic.h>
 *
 * @brief Base concept for streams.
 *
 * @signature concept StreamConcept;
 */

/*!
 * @mfn StreamConcept#Value
 * @brief Metafunction for retrieving the value type of a stream.
 *
 * @signature Value<TStream>::Type;
 *
 * @tparam TStream The stream type to query for its value type.
 * @return Type    The resulting value type.
 */

/*!
 * @mfn StreamConcept#Size
 * @brief Metafunction for retrieving the type of a stream.
 *
 * @signature Size<TStream>::Type;
 *
 * @tparam TStream The stream type to query for its size type.
 * @return Type    The resulting size type.
 */

/*!
 * @mfn StreamConcept#Position
 * @brief Metafunction for retrieving the position type of a stream.
 *
 * @signature Position<TStream>::Type;
 *
 * @tparam TStream The stream type to query for its position type.
 * @return Type    The resulting position type.
 */

/*!
 * @fn StreamConcept#position
 * @brief Return current stream position.
 *
 * @signature TPosition position(stream);
 *
 * @param[in] stream The stream to query.
 * @return TPosition Current position in stream, see @link StreamConcept#Position Position @endlink.
 */

/*!
 * @fn StreamConcept#setPosition
 * @brief Set stream position.
 *
 * @signature void setPosition(stream, pos);
 *
 * @param[in,out] stream The stream to update
 * @param[in]     pos    The positoin to set.
 */

/*!
 * @fn StreamConcept#atEnd
 * @brief Return whether stream is at the end.
 *
 * @signature bool atEnd(stream);
 *
 * @param[in] stream The stream to check.
 * @return bool <tt>true</tt> if the file at EOF, <tt>false</tt> otherwise.
 */

SEQAN_CONCEPT(StreamConcept, (TStream))
{};

// --------------------------------------------------------------------------
// Concept InputStreamConcept
// --------------------------------------------------------------------------

/*!
 * @concept InputStreamConcept Input StreamConcept
 * @extends StreamConcept
 * @headerfile <seqan/basic.h>
 *
 * @signature concept InputStreamConcept : StreamConcept;
 *
 * @brief Concept for input streams (for reading).
 */

SEQAN_CONCEPT_REFINE(InputStreamConcept, (TStream), (StreamConcept))
{
    typedef typename Value<TStream>::Type       TValue;
    typedef typename Size<TStream>::Type        TSize;
    typedef typename Position<TStream>::Type    TPosition;

    SEQAN_CONCEPT_ASSERT((SignedIntegerConcept<TPosition>));

    SEQAN_CONCEPT_USAGE(InputStreamConcept)
    {}
};

// --------------------------------------------------------------------------
// Concept OutputStreamConcept
// --------------------------------------------------------------------------

/*!
 * @concept OutputStreamConcept Output StreamConcept
 * @extends StreamConcept
 * @headerfile <seqan/basic.h>
 *
 * @signature concept OutputStreamConcept : StreamConcept;
 *
 * @brief Concept for output streams (for writing).
 */

SEQAN_CONCEPT_REFINE(OutputStreamConcept, (TStream), (StreamConcept))
{
    typedef typename Value<TStream>::Type       TValue;
    typedef typename Size<TStream>::Type        TSize;
    typedef typename Position<TStream>::Type    TPosition;

    SEQAN_CONCEPT_ASSERT((SignedIntegerConcept<TPosition>));

    SEQAN_CONCEPT_USAGE(OutputStreamConcept)
    {}
};

// --------------------------------------------------------------------------
// Concept BidirectionalStreamConcept
// --------------------------------------------------------------------------

/*!
 * @concept BidirectionalStreamConcept Bidirectional StreamConcept
 * @extends StreamConcept
 * @headerfile <seqan/basic.h>
 *
 * @signature concept BidirectionalStreamConcept : StreamConcept;
 *
 * @brief Concept for bidirectional streams (both for reading and writing).
 */

SEQAN_CONCEPT_REFINE(BidirectionalStreamConcept, (TStream), (InputStreamConcept)(OutputStreamConcept))
{};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Struct FormattedNumber
// ----------------------------------------------------------------------------

/*!
 * @class FormattedNumber
 * @implements NumberConcept
 * @headerfile <seqan/basic.h>
 *
 * @brief Helper class for storing a numeric value together with a
 *        <a href="http://www.cplusplus.com/reference/cstdio/printf/">printf format string</a>.
 *
 * @signature template <typename TValue>
 *            struct FormattedNumber;
 *
 * @tparam The numeric value type.
 */

template <typename TValue>
struct FormattedNumber
{
    char const * format;
    TValue value;

    /*!
     * @fn FormattedNumber::FormattedNumber
     * @brief Constructor.
     *
     * @signature FormattedNumber::FormattedNumber(format, value);
     *
     * @param[in] format A <tt>char const *</tt> for the format string.
     * @param[in] value  The <tt>TValue</tt> to store.
     *
     * The constructed FormattedNumber object store the <tt>format</tt> pointer "as is".  This means that you are
     * responsible for keeping this pointer valid until the object is deconstructed.  Passing in a C string literal
     * (as in <tt>FormattedNumber&lt;double&gt;("%.2f", 1.234)</tt> is fine.
     */

    FormattedNumber(char const * format, TValue const & value) :
        format(format), value(value)
    {}

    operator TValue() const
    {
        return value;
    }
};

template <typename TValue>
struct Is< NumberConcept< FormattedNumber<TValue> > > :
    Is< NumberConcept<TValue> > {};

// ============================================================================
// Exceptions
// ============================================================================

// ----------------------------------------------------------------------------
// Exception ParseError
// ----------------------------------------------------------------------------

/*!
 * @class ParseError
 * @headerfile <seqan/basic.h>
 *
 * @brief Exception class for parser errors.
 *
 * @signature struct ParserError : RuntimeError;
 */

struct ParseError : RuntimeError
{
    /*!
     * @fn ParseError::ParseError
     * @headerfile <seqan/basic.h>
     *
     * @brief Constructor.
     *
     * @signature ParseError::ParseError(message);
     *
     * @param[in] message The error message to use, <tt>std::string</tt> or <tt>char const * </tt>.
     */

    template <typename TString>
    ParseError(TString const & message) :
        RuntimeError(message)
    {}
};

// ----------------------------------------------------------------------------
// Exception UnexpectedEnd
// ----------------------------------------------------------------------------

/*!
 * @class UnexpectedEnd
 * @headerfile <seqan/basic.h>
 *
 * @brief Exception class for "unexpected end of input" errors.
 *
 * @signature struct UnexpectedEnd : RuntimeError;
 */

struct UnexpectedEnd : ParseError
{
    /*!
     * @fn UnexpectedEnd::UnexpectedEnd
     * @headerfile <seqan/basic.h>
     *
     * @brief Default constructor, makes the object use a default message.
     *
     * @signature UnexpectedEnd::UnexpectedEnd();
     */

    UnexpectedEnd() :
        ParseError("Unexpected end of input.")
    {}
};

// ----------------------------------------------------------------------------
// Exception EmptyFieldError
// ----------------------------------------------------------------------------

/*!
 * @class EmptyFieldError
 * @headerfile <seqan/basic.h>
 *
 * @brief Exception class for "empty field" errors.
 *
 * @signature struct EmptyFieldError : RuntimeError;
 */

struct EmptyFieldError : ParseError
{
    /*!
     * @fn EmptyFieldError::EmptyFieldError
     * @headerfile <seqan/basic.h>
     *
     * @brief Construct the exception with <tt>fieldName + " field was empty."</tt>.
     *
     * @signature EmptyFieldEror::EmptyFieldError(fieldName);
     *
     * @param[in] fieldName The field name to use for the message, <tt>std::string</tt>.
     */

    EmptyFieldError(std::string fieldName):
        ParseError(fieldName + " field was empty.")
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

/*!
 * @mfn StreamConcept#DirectionIterator
 * @brief Return the direction iterator for the given direction.
 *
 * @signature DirectionIterator<TStream>::Type;
 *
 * @tparam TStream The stream to query for its direction iterator.
 * @return Type    The resulting direction iterator.
 */

/*!
 * @mfn ContainerConcept#DirectionIterator
 * @brief Return the direction iterator for the given direction.
 *
 * @signature DirectionIterator<TContainer>::Type;
 *
 * @tparam TContainer The container to query for its direction iterator.
 * @return Type       The resulting direction iterator.
 */

template <typename TObject, typename TDirection>
struct DirectionIterator :
    If<Is<StreamConcept<TObject> >,
       Iter<TObject, StreamIterator<TDirection> >,
       typename Iterator<TObject, Rooted>::Type>
{};

// --------------------------------------------------------------------------
// Metafunction BasicStream
// --------------------------------------------------------------------------

/*!
 * @mfn BasicStream
 * @headerfile <seqan/basic.h>
 * @brief Return the stream type to read or write values.
 *
 * @signature BasicStream<TValue, TDirection[, TTraits]>::Type;
 *
 * @tparam TValue     The value type of the stream.
 * @tparam TDirection The direction of the stream, one of the @link DirectionTags @endlink.
 * @tparam TTraits    The traits to use for the values, defaults to <tt>std::char_traits&lt;TValue&gt;</tt>.
 *
 */

template <typename TValue, typename TDirection, typename TTraits = std::char_traits<TValue> >
struct BasicStream :
    If<
        IsSameType<TDirection, Input>,
        std::basic_istream<TValue, TTraits>,
        typename If<
            IsSameType<TDirection, Output>,
            std::basic_ostream<TValue, TTraits>,
            std::basic_iostream<TValue, TTraits>
            >::Type
        >
{};

// --------------------------------------------------------------------------
// Metafunction IosOpenMode
// --------------------------------------------------------------------------

/*!
 * @mfn IosOpenMode
 * @headerfile <seqan/basic.h>
 * @brief Return the <tt>std::ios</tt> open mode for a direction.
 *
 * @signature IosOpenMode<TDirection[, TDummy]>::Type;
 *
 * @tparam TDirection The direction to query for the open mode, one of the @link DirectionTags @endlink.
 * @tparam TDummy     Implementation detail, defaults to <tt>void</tt> and is ignored.
 * @return Type       The resulting open mode of type <tt>const int</tt>.
 */

template <typename TDirection, typename TDummy = void>
struct IosOpenMode;


template <typename TDummy>
struct IosOpenMode<Input, TDummy>
{
    static const int VALUE;
};

template <typename TDummy>
struct IosOpenMode<Output, TDummy>
{
    static const int VALUE;
};

template <typename TDummy>
struct IosOpenMode<Bidirectional, TDummy>
{
    static const int VALUE;
};

template <typename TDummy>
const int IosOpenMode<Input, TDummy>::VALUE = std::ios::in | std::ios::binary;

template <typename TDummy>
const int IosOpenMode<Output, TDummy>::VALUE = std::ios::out | std::ios::binary;

template <typename TDummy>
const int IosOpenMode<Bidirectional, TDummy>::VALUE = std::ios::in | std::ios::out | std::ios::binary;


// --------------------------------------------------------------------------
// Metafunction MagicHeader
// --------------------------------------------------------------------------

/*!
 * @mfn MagicHeader
 * @headerfile <seqan/basic.h>
 * @brief Returns the magic header for a file format tag.
 *
 * The magic header is used for recognizing files from the first few bytes.
 *
 * @signature MagicHeader<TTag[, TDummy]>::VALUE;
 *
 * @tparam TTag   The file format tag to use for the query.
 * @tparam TDummy Implementation detail, defaults to <tt>void</tt> and is ignored.
 * @return VALUE  The magic header string, of type <tt>char const *</tt>.
 *
 * This metafunction must be implemented in the modules implementing the file I/O.  The metafunction is predefined when
 * <tt>TTag</tt> is @link Nothing @endlink.  In this case, <tt>VALUE</tt> is <tt>NULL</tt>.
 */

template <typename TTag, typename T = void>
struct MagicHeader;

template <typename T>
struct MagicHeader<Nothing, T>
{
    static char const * VALUE;
};

template <typename T>
char const * MagicHeader<Nothing, T>::VALUE = NULL;

// --------------------------------------------------------------------------
// Metafunction FileExtensions
// --------------------------------------------------------------------------

/*!
 * @mfn FileExtensions
 * @headerfile <seqan/basic.h>
 * @brief Returns an array of file format extension strings for file foramt tag.
 *
 * @signature FileExtensions<TFormat[, TDummy]>::VALUE;
 *
 * @tparam TTag   The file format tag to use for the query.
 * @tparam TDummy Implementation detail, defaults to <tt>void</tt> and is ignored.
 * @return VALUE  The array of file format extension, of type <tt>char const *[]</tt>.
 *
 * This metafunction must be implemented in the modules implementing the file I/O.  The metafunction is predefined when
 * <tt>TTag</tt> is @link Nothing @endlink.  In this case, <tt>VALUE</tt> is <tt>{""}</tt>.
 */

template <typename TFormat, typename T = void>
struct FileExtensions;

template <typename T>
struct FileExtensions<Nothing, T>
{
    static char const * VALUE[1];
};

template <typename T>
char const * FileExtensions<Nothing, T>::VALUE[1] =
{
    ""  // default output extension
};

// ----------------------------------------------------------------------------
// Metafunction IntegerFormatString_
// ----------------------------------------------------------------------------
// Return the format string for numbers.

template <typename TUnsigned, unsigned SIZE, typename T = void>
struct IntegerFormatString_;


template <typename TUnsigned, typename T>
struct IntegerFormatString_<TUnsigned, 1, T> :
    IntegerFormatString_<TUnsigned, 2, T> {};


template <typename T>
struct IntegerFormatString_<False, 2, T>
{
    static const char VALUE[];
    typedef short Type;
};
template <typename T>
const char IntegerFormatString_<False, 2, T>::VALUE[] = "%hi";


template <typename T>
struct IntegerFormatString_<True, 2, T>
{
    static const char VALUE[];
    typedef unsigned short Type;
};
template <typename T>
const char IntegerFormatString_<True, 2, T>::VALUE[] = "%hu";


template <typename T>
struct IntegerFormatString_<False, 4, T>
{
    static const char VALUE[];
    typedef int Type;
};
template <typename T>
const char IntegerFormatString_<False, 4, T>::VALUE[] = "%i";


template <typename T>
struct IntegerFormatString_<True, 4, T>
{
    static const char VALUE[];
    typedef unsigned Type;
};
template <typename T>
const char IntegerFormatString_<True, 4, T>::VALUE[] = "%u";


// helper for the case: typedef long __int64;
template <typename TIsUnsigned, typename T>
struct LongFormatString_;

template <typename T>
struct LongFormatString_<False, T>
{
    static const char VALUE[];
    typedef long Type;
};
template <typename T>
const char LongFormatString_<False, T>::VALUE[] = "%li";

template <typename T>
struct LongFormatString_<True, T>
{
    static const char VALUE[];
    typedef unsigned long Type;
};
template <typename T>
const char LongFormatString_<True, T>::VALUE[] = "%lu";

// helper for the case: typedef long long __int64;
template <typename TIsUnsigned, typename T>
struct Int64FormatString_;

template <typename T>
struct Int64FormatString_<False, T>
{
    static const char VALUE[];
    typedef __int64 Type;
};
template <typename T>
const char Int64FormatString_<False, T>::VALUE[] = "%lli";

template <typename T>
struct Int64FormatString_<True, T>
{
    static const char VALUE[];
    typedef __uint64 Type;
};
template <typename T>
const char Int64FormatString_<True, T>::VALUE[] = "%llu";


template <typename TIsUnsigned, typename T>
struct IntegerFormatString_<TIsUnsigned, 8, T> :
    If<IsSameType<__uint64, unsigned long>,
       LongFormatString_<TIsUnsigned, T>,
       Int64FormatString_<TIsUnsigned, T> >::Type {};


// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function writeValue()                                     [ContainerConcept]
// ----------------------------------------------------------------------------

/*!
 * @fn ContainerConcept#writeValue
 * @brief Write a value at the end of a container.
 *
 * @signature void writeValue(container, val);
 *
 * @param[in,out] container to append to.
 * @param[in]     val       The value to append.
 *
 * @see ContainerConcept#appendValue
 */

// resizable containers
template <typename TSequence, typename TValue>
inline SEQAN_FUNC_ENABLE_IF(Is<ContainerConcept<TSequence> >, void)
writeValue(TSequence &cont, TValue val)
{
    appendValue(cont, val);
}

// ----------------------------------------------------------------------------
// Function writeValue()                                                [Range]
// ----------------------------------------------------------------------------

/*!
 * @fn Range#writeValue
 * @brief Write a value to a @link Range @endlink.
 *
 * @signature void writeValue(range, val);
 *
 * <tt>val</tt> will be assigned to the first element of the range.  Then, the beginning of the range will be advanced
 * by one.
 *
 * @param[in,out] range to append to.
 * @param[in]     val   The value to append.
 */

// Range
template <typename TIterator, typename TValue>
inline void
writeValue(Range<TIterator> &range, TValue val)
{
    assignValue(range.begin, val);
    ++range.begin;
}

// ----------------------------------------------------------------------------
// Function writeValue()                                                 [Iter]
// ----------------------------------------------------------------------------

/*!
 * @fn OutputIteratorConcept#writeValue
 * @brief Write a single value to a container by dereferencing its iterator.
 *
 * @signature void writeValue(iter, val);
 *
 * @param[in,out] iter The iterator to use for dereferenced writing.
 * @param[in]     val  The value to write into the container.
 *
 * If the host of <tt>iter</tt> is a @link ContainerConcept @endlink then container is resized to make space for the
 * item.
 */

// resizable containers
template <typename TSequence, typename TSpec, typename TValue>
inline SEQAN_FUNC_ENABLE_IF(Is<ContainerConcept<TSequence> >, void)
writeValue(Iter<TSequence, TSpec> & iter, TValue val)
{
    typedef Iter<TSequence, TSpec> TIter;

    TSequence &cont = container(iter);
    typename Position<TIter>::Type pos = position(iter);
    typename Size<TIter>::Type len = length(cont);

    if (pos < len)
    {
        assignValue(iter, val);
        ++iter;
    }
    else
    {
        if (pos > len)
            resize(cont, pos - 1);
        appendValue(cont, val);
        setPosition(iter, pos + 1);
    }
}

// non-resizable containers
template <typename TNoSequence, typename TSpec, typename TValue>
inline SEQAN_FUNC_DISABLE_IF(Is<ContainerConcept<TNoSequence> >, void)
writeValue(Iter<TNoSequence, TSpec> & iter, TValue val)
{
    SEQAN_ASSERT_LT(position(iter), length(container(iter)));

    assignValue(iter, val);
    ++iter;
}

// ----------------------------------------------------------------------------
// Function writeValue()                                              [pointer]
// ----------------------------------------------------------------------------

///*!
// * @fn ContainerConcept#writeValue
// * @brief Write a value by dereferencing a pointer and incrementing its position by one.
// *
// * @signature void writeValue(pointer, val);
// *
// * @param[in,out] iter The pointer to dereference, usually a <tt>char *</tt>.
// * @param[in]     val  The value to write to the dereferenced pointer.
// *
// * This function is equivalent to <tt>*iter++ = val</tt>.
// */

template <typename TTargetValue, typename TValue>
inline void
writeValue(TTargetValue * & iter, TValue val)
{
    *iter++ = val;
}

// ----------------------------------------------------------------------------
// Function _write(); Element-wise
// ----------------------------------------------------------------------------

template <typename TTarget, typename TFwdIterator, typename TSize, typename TIChunk, typename TOChunk>
inline void _write(TTarget &target, TFwdIterator &iter, TSize n, TIChunk, TOChunk)
{
    for (; n > (TSize)0; --n, ++iter)
        writeValue(target, getValue(iter));
}

// ----------------------------------------------------------------------------
// Function _write(); Chunked
// ----------------------------------------------------------------------------

template <typename TTarget, typename TFwdIterator, typename TSize, typename TIValue, typename TOValue>
inline void _write(TTarget &target, TFwdIterator &iter, TSize n, Range<TIValue*> *, Range<TOValue*> *)
{
    typedef Nothing* TNoChunking;
    typedef typename Size<TTarget>::Type TTargetSize;

    Range<TIValue*> ichunk;
    Range<TOValue*> ochunk;

    while (n != (TSize)0)
    {
        getChunk(ichunk, iter, Input());
        getChunk(ochunk, target, Output());

        TTargetSize minChunkSize = std::min((TTargetSize)length(ichunk), (TTargetSize)length(ochunk));

        if (SEQAN_UNLIKELY(minChunkSize == 0u))
        {
            reserveChunk(target, n, Output());
            reserveChunk(iter, n, Input());
            getChunk(ochunk, target, Output());
            getChunk(ichunk, iter, Input());
            minChunkSize = std::min((TTargetSize)length(ichunk), (TTargetSize)length(ochunk));
            if (SEQAN_UNLIKELY(minChunkSize == 0u))
            {
                _write(target, iter, n, TNoChunking(), TNoChunking());
                return;
            }
        }

        if (minChunkSize > (TTargetSize)n)
            minChunkSize = (TTargetSize)n;

        arrayCopyForward(ichunk.begin, ichunk.begin + minChunkSize, ochunk.begin);

        iter += minChunkSize;                      // advance input iterator
        advanceChunk(target, minChunkSize);
        n -= minChunkSize;
    }
}

// chunked, target is pointer (e.g. readRawPod)
template <typename TOValue, typename TFwdIterator, typename TSize>
inline SEQAN_FUNC_DISABLE_IF(IsSameType<typename Chunk<TFwdIterator>::Type, Nothing>, void)
write(TOValue *ptr, TFwdIterator &iter, TSize n)
{
    typedef Nothing* TNoChunking;
    typedef typename Size<TFwdIterator>::Type TSourceSize;
    typedef typename Chunk<TFwdIterator>::Type TIChunk;

    TIChunk ichunk;

    while (n != (TSize)0)
    {
        getChunk(ichunk, iter, Input());
        TSourceSize chunkSize = length(ichunk);

        if (SEQAN_UNLIKELY(chunkSize == 0u))
        {
            reserveChunk(iter, n, Input());
            getChunk(ichunk, iter, Input());
            TSourceSize chunkSize = length(ichunk);
            if (SEQAN_UNLIKELY(chunkSize == 0u))
            {
                _write(ptr, iter, n, TNoChunking(), TNoChunking());
                return;
            }
        }

        if (chunkSize > (TSourceSize)n)
            chunkSize = (TSourceSize)n;

        arrayCopyForward(ichunk.begin, ichunk.begin + chunkSize, ptr);

        iter += chunkSize;                          // advance input iterator
        ptr += chunkSize;
        n -= chunkSize;
    }
}

// non-chunked fallback
template <typename TTarget, typename TIValue, typename TSize>
inline SEQAN_FUNC_ENABLE_IF(And< IsSameType<typename Chunk<TTarget>::Type, Nothing>,
                                 Is<Convertible<typename Value<TTarget>::Type, TIValue> > >, void)
write(TTarget &target, TIValue *ptr, TSize n)
{
    _write(target, ptr, n, Nothing(), Nothing());
}

// ostream shortcut, source is pointer (e.g. readRawPod)
template <typename TTarget, typename TSize>
inline SEQAN_FUNC_ENABLE_IF(Is< OutputStreamConcept<TTarget> >, void)
write(TTarget &target, const char *ptr, TSize n)
{
    target.write(ptr, n);
}

// ostream shortcut, source is pointer (e.g. readRawPod)
template <typename TTarget, typename TSize>
inline SEQAN_FUNC_ENABLE_IF(Is< OutputStreamConcept<TTarget> >, void)
write(TTarget &target, char *ptr, TSize n)
{
    target.write(ptr, n);
}

// chunked, source is pointer (e.g. readRawPod)
template <typename TTarget, typename TIValue, typename TSize>
inline SEQAN_FUNC_ENABLE_IF(And< Not<IsSameType<typename Chunk<TTarget>::Type, Nothing> >,
                                 Is<Convertible<typename Value<TTarget>::Type, TIValue> > >, void)
write(TTarget &target, TIValue *ptr, TSize n)
{
    typedef Nothing* TNoChunking;
    typedef typename Size<TTarget>::Type TTargetSize;
    typedef typename Chunk<TTarget>::Type TOChunk;

    TOChunk ochunk;

    while (n != (TSize)0)
    {
        getChunk(ochunk, target, Output());
        TTargetSize chunkSize = length(ochunk);

        if (SEQAN_UNLIKELY(chunkSize == 0u))
        {
            reserveChunk(target, n, Output());
            getChunk(ochunk, target, Output());
            chunkSize = length(ochunk);
            if (SEQAN_UNLIKELY(chunkSize == 0u))
            {
                _write(target, ptr, n, TNoChunking(), TNoChunking());
                return;
            }
        }

        if (chunkSize > (TTargetSize)n)
            chunkSize = (TTargetSize)n;

        arrayCopyForward(ptr, ptr + chunkSize, ochunk.begin);

        ptr += chunkSize;                      // advance input iterator
        advanceChunk(target, chunkSize);
        n -= chunkSize;
    }
}

template <typename TOValue, typename TIValue, typename TSize>
inline SEQAN_FUNC_ENABLE_IF(And< Is<CharConcept<TOValue> >,
                                 Is<CharConcept<TIValue> > >, void)
write(TOValue * &optr, TIValue *iptr, TSize n)
{
    std::memcpy(optr, iptr, n);
    optr += n;
}

template <typename TOValue, typename TIValue, typename TSize>
inline SEQAN_FUNC_ENABLE_IF(And< Is<CharConcept<TOValue> >,
                                 Is<CharConcept<TIValue> > >, void)
write(TOValue * optr, TIValue * &iptr, TSize n)
{
    std::memcpy(optr, iptr, n);
    iptr += n;
}

// ----------------------------------------------------------------------------
// Function write(TValue *)
// ----------------------------------------------------------------------------
// NOTE(esiragusa): should it be defined for Streams and Containers?

//template <typename TTarget, typename TValue, typename TSize>
//inline SEQAN_FUNC_ENABLE_IF(Or<Is<OutputStreamConcept<TTarget> >, Is<ContainerConcept<TTarget> > >, void)
//write(TTarget &target, TValue *ptr, TSize n)
//{
//    typedef Range<TValue*>                          TRange;
//    typedef typename Iterator<TRange, Rooted>::Type TIterator;
//    typedef typename Chunk<TIterator>::Type*        TIChunk;
//    typedef typename Chunk<TTarget>::Type*          TOChunk;
//
//    TRange range(ptr, ptr + n);
//    TIterator iter = begin(range, Rooted());
//    _write(target, iter, n, TIChunk(), TOChunk());
//}

// ----------------------------------------------------------------------------
// Function write(Iterator<Input>)
// ----------------------------------------------------------------------------

/*!
 * @fn ContainerConcept#write
 * @brief Write to a container.
 *
 * @signature void write(container, iter, n);
 *
 * @param[in,out] container The container to append to.
 * @param[in,out] iter      The @link ForwardIteratorConcept forward iterator @endlink to take the values from.
 * @param[in]     n         Number of elements to write from <tt>iter</tt>.
 *
 * This function reads <tt>n</tt> values from <tt>iter</tt> and appends them to the back of <tt>container</tt>.
 */

//TODO(singer): Enable this!
template <typename TTarget, typename TFwdIterator, typename TSize>
//inline SEQAN_FUNC_ENABLE_IF(Or<Is<OutputStreamConcept<TTarget> >, Is<ContainerConcept<TTarget> > >, void)
inline SEQAN_FUNC_ENABLE_IF(And< Is<IntegerConcept<TSize> >,
                                 Is<Convertible<typename Value<TTarget>::Type,
                                                typename Value<TFwdIterator>::Type> > >, void)
write(TTarget &target, TFwdIterator &iter, TSize n)
{
    typedef typename Chunk<TFwdIterator>::Type* TIChunk;
    typedef typename Chunk<TTarget>::Type*      TOChunk;

    _write(target, iter, n, TIChunk(), TOChunk());
}

// write for more complex values (defer to write of iterator value)
// used for Strings of Pairs
template <typename TTarget, typename TFwdIterator, typename TSize>
//inline SEQAN_FUNC_ENABLE_IF(Or<Is<OutputStreamConcept<TTarget> >, Is<ContainerConcept<TTarget> > >, void)
inline SEQAN_FUNC_ENABLE_IF(And<
    Is<IntegerConcept<TSize> >,
    Not< Is<Convertible<typename Value<TTarget>::Type,
                        typename Value<TFwdIterator>::Type> > > >, void)
write(TTarget &target, TFwdIterator &iter, TSize n)
{
    for (; n > (TSize)0; --n, ++iter)
    {
        write(target, *iter);
        writeValue(target, ' ');
    }
}

// ----------------------------------------------------------------------------
// Function write(TContainer)
// ----------------------------------------------------------------------------

template <typename TTarget, typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(And< Is<ContainerConcept<TContainer> >,
                                 Not<IsContiguous<TContainer> > >, void)
write(TTarget &target, TContainer &cont)
{
    typename DirectionIterator<TContainer, Input>::Type iter = directionIterator(cont, Input());
    write(target, iter, length(cont));
}

template <typename TTarget, typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(And< Is<ContainerConcept<TContainer> >,
                                 IsContiguous<TContainer> >, void)
write(TTarget &target, TContainer &cont)
{
    typename Iterator<TContainer, Standard>::Type iter = begin(cont, Standard());
    write(target, iter, length(cont));
}



template <typename TTarget, typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(And< Is<ContainerConcept<TContainer> >,
                                 Not<IsContiguous<TContainer> > >, void)
write(TTarget &target, TContainer const &cont)
{
    typename DirectionIterator<TContainer const, Input>::Type iter = directionIterator(cont, Input());
    write(target, iter, length(cont));
}

template <typename TTarget, typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(And< Is<ContainerConcept<TContainer> >,
                                 IsContiguous<TContainer> >, void)
write(TTarget &target, TContainer const &cont)
{
    typename Iterator<TContainer const, Standard>::Type iter = begin(cont, Standard());
    write(target, iter, length(cont));
}



template <typename TTarget, typename TValue>
inline void
write(TTarget &target, TValue * ptr)
{
    write(target, ptr, length(ptr));
}

// ----------------------------------------------------------------------------
// Function appendNumber()
// ----------------------------------------------------------------------------
// Generic version for integers.

template <typename TTarget, typename TInteger>
inline SEQAN_FUNC_ENABLE_IF(Is<IntegerConcept<TInteger> >, typename Size<TTarget>::Type)
appendNumber(TTarget & target, TInteger i)
{
    typedef IntegerFormatString_<typename Is<UnsignedIntegerConcept<TInteger> >::Type,
                                 sizeof(TInteger)> TInt;

    // 1 byte has at most 3 decimal digits (plus 2 for '-' and the NULL character)
    char buffer[sizeof(TInteger) * 3 + 2];
    size_t len = snprintf(buffer, sizeof(buffer),
                          TInt::VALUE, static_cast<typename TInt::Type>(i));
    char *bufPtr = buffer;
    write(target, bufPtr, len);
    return len;
}

// ----------------------------------------------------------------------------
// Function appendNumber(bool)
// ----------------------------------------------------------------------------

template <typename TTarget>
inline typename Size<TTarget>::Type
appendNumber(TTarget & target, bool source)
{
    writeValue(target, '0' + source);
    return 1;
}

// ----------------------------------------------------------------------------
// Function appendNumber(float)
// ----------------------------------------------------------------------------

template <typename TTarget>
inline typename Size<TTarget>::Type
appendNumber(TTarget & target, float source)
{
    char buffer[32];
    size_t len = snprintf(buffer, sizeof(buffer), "%g", source);
    write(target, (char *)buffer, len);
    return len;
}

// ----------------------------------------------------------------------------
// Function appendNumber(double)
// ----------------------------------------------------------------------------

template <typename TTarget>
inline typename Size<TTarget>::Type
appendNumber(TTarget & target, double source)
{
    char buffer[32];
    size_t len = snprintf(buffer, sizeof(buffer), "%g", source);
    write(target, (char *)buffer, len);
    return len;
}

// ----------------------------------------------------------------------------
// Function appendNumber(double)
// ----------------------------------------------------------------------------

template <typename TTarget, typename TValue>
inline typename Size<TTarget>::Type
appendNumber(TTarget & target, FormattedNumber<TValue> const & source)
{
    char buffer[100];
    size_t len = snprintf(buffer, sizeof(buffer), source.format, source.value);
    write(target, (char *)buffer, len);
    return len;
}

template <typename TValue>
inline FormattedNumber<TValue>
formattedNumber(const char *format, TValue const & val)
{
    return FormattedNumber<TValue>(format, val);
}

// ----------------------------------------------------------------------------
// Function appendRawPod()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TValue>
inline void
appendRawPod(TTarget & target, TValue const & val)
{
    write(target, (unsigned char*)&val, sizeof(TValue));
}

template <typename TTargetValue, typename TValue>
inline void
appendRawPod(TTargetValue * &ptr, TValue const & val)
{
    *reinterpret_cast<TValue* &>(ptr)++ = val;
}

// ----------------------------------------------------------------------------
// Function write(TNumber); write fundamental type
// ----------------------------------------------------------------------------

template <typename TTarget, typename TValue>
inline SEQAN_FUNC_ENABLE_IF(And< Is<Convertible<typename Value<TTarget>::Type, TValue> >,
                                 Is<FundamentalConcept<TValue> > >, void)
write(TTarget &target, TValue &number)
{
    if (sizeof(TValue) == 1)
        writeValue(target, number);     // write chars as chars
    else
        appendNumber(target, number);
}

template <typename TTarget, typename TValue>
inline SEQAN_FUNC_ENABLE_IF(And< Is<Convertible<typename Value<TTarget>::Type, TValue> >,
                                 Is<FundamentalConcept<TValue> > >, void)
write(TTarget &target, TValue const &number)
{
    if (sizeof(TValue) == 1)
        writeValue(target, number);     // write chars as chars
    else
        appendNumber(target, number);
}

// ----------------------------------------------------------------------------
// Function write(TNumber); write non-fundamental, convertible type
// ----------------------------------------------------------------------------

template <typename TTarget, typename TValue>
inline SEQAN_FUNC_ENABLE_IF(And< Is<Convertible<typename Value<TTarget>::Type, TValue> >,
                                 Not<Is<FundamentalConcept<TValue> > > >, void)
write(TTarget &target, TValue &number)
{
    writeValue(target, number);
}

template <typename TTarget, typename TValue>
inline SEQAN_FUNC_ENABLE_IF(And< Is<Convertible<typename Value<TTarget>::Type, TValue const> >,
                                 Not<Is<FundamentalConcept<TValue const> > > >, void)
write(TTarget &target, TValue const &number)
{
    writeValue(target, number);
}

// ----------------------------------------------------------------------------
// Function read(Iterator<Input>)
// ----------------------------------------------------------------------------

template <typename TTarget, typename TFwdIterator, typename TSize>
inline SEQAN_FUNC_ENABLE_IF(Is<IntegerConcept<TSize> >, TSize)
read(TTarget &target, TFwdIterator &iter, TSize n)
{
    TSize i;
    for (i = 0; !atEnd(iter) && i < n; ++i, ++iter)
        writeValue(target, value(iter));
    return i;
}

// ----------------------------------------------------------------------------
// Function read(TContainer)
// ----------------------------------------------------------------------------

template <typename TTarget, typename TContainer>
inline typename Size<TTarget>::Type
read(TTarget &target, TContainer &cont)
{
    typename DirectionIterator<TContainer, Input>::Type iter = directionIterator(cont, Input());
    return read(target, iter, length(cont));
}

// ----------------------------------------------------------------------------
// operator<<
// ----------------------------------------------------------------------------

template <typename TStream, typename TIterator>
inline TStream &
operator<<(TStream & target,
           Range<TIterator> const & source)
{
    typename DirectionIterator<TStream, Output>::Type it = directionIterator(target, Output());
    write(it, source);
    return target;
}

template <typename TStream, typename TValue>
inline TStream &
operator<<(TStream & target,
           FormattedNumber<TValue> const & source)
{
    typename DirectionIterator<TStream, Output>::Type it = directionIterator(target, Output());
    write(it, source);
    return target;
}

}  // namespace seqean

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_BASIC_STREAM_H_
