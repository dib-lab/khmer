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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// String => Numerical conversions
// ==========================================================================

#ifndef SEQAN_STREAM_LEXICAL_CAST_H
#define SEQAN_STREAM_LEXICAL_CAST_H

namespace seqan {

// ============================================================================
// Exceptions
// ============================================================================

// ----------------------------------------------------------------------------
// Exception BadLexicalCast
// ----------------------------------------------------------------------------

/*!
 * @class BadLexicalCast
 * @extends ParseError
 * @headerfile <seqan/stream.h>
 * @brief Throw on bad lexical casts.
 *
 * @signature struct BadLexicalCast : ParseError;
 */

struct BadLexicalCast : ParseError
{
    /*!
     * @fn BadLexicalCast::BadLexicalCast
     * @brief Constructor.
     *
     * @signature BadLexicalCast::BadLexicalCast(target, source);
     *
     * @param[in] target Target value, used as a tag only.
     * @param[in] source Source value, a @link StringConcept sequence @endlink of <tt>char</tt>.
     */
    template <typename TTarget, typename TSource>
    BadLexicalCast(TTarget const & target, TSource const & source) :
        ParseError(std::string("Unable to convert '") +
                   std::string(begin(source, Standard()), end(source, Standard())) +
                   "' into " + toCString(Demangler<TTarget>(target)) + ".")
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function lexicalCast()
// ----------------------------------------------------------------------------

/*!
 * @fn lexicalCast
 * @headerfile <seqan/stream.h>
 * @brief Interpret the character sequence in <tt>source</tt> as the numeric value and write to <tt>target</tt>.
 *
 * @signature bool lexicalCast(target, source);
 * @signature TTarget lexicalCast<TTarge>(source)  // throws BadLexicalCast
 *
 * @param[out] target  A numeric value to write to.
 * @param[in]  source  The source sequence of <tt>char</tt> to convert.
 * @tparam     TTarget The type to use for lexical cast of <tt>source</tt>.
 *
 * @throw BadLexicalCast The second variant throws @link BadLexicalCast @endlink in the case that casting failed.
 *
 * @return bool <tt>true</tt> if successful and <tt>false</tt> if there was a problem with the cast.
 *
 * @section Examples
 *
 * The following example shows some lexical cast from various sequence types to numbers.
 *
 * @include demos/stream/lexical_cast.cpp
 *
 * @include demos/stream/lexical_cast.cpp.stdout
 */

// single char
template <typename TSource>
inline bool
lexicalCast(char & target, TSource const & source)
{
    if (SEQAN_UNLIKELY(length(source) != 1))
        return false;
    target = getValue(begin(source, Standard()));
    return true;
}

// (weese:) we have to implement our own cast functions as not all sources support toCString()

// Generic version for unsigned integers.
template <typename TInteger, typename TSource>
inline SEQAN_FUNC_ENABLE_IF(Is<UnsignedIntegerConcept<TInteger> >, bool)
lexicalCast(TInteger & target, TSource const & source)
{
    typedef typename Iterator<TSource const, Standard>::Type TIter;

    TIter it = begin(source, Standard());
    TIter itEnd = end(source, Standard());

    if (SEQAN_UNLIKELY(it == itEnd))
        return false;

    TInteger val = 0;
    do
    {
        unsigned char digit = *it++ - '0';

        // invalid digit detection
        if (SEQAN_UNLIKELY(digit > 9))
            return false;

        // overflow detection
        if (SEQAN_UNLIKELY(val > MaxValue<TInteger>::VALUE / 10))
            return false;
        val *= 10;

        // overflow detection
        val += digit;
        if (SEQAN_UNLIKELY(val < digit))
            return false;
    }
    while (it != itEnd);
    target = val;
    return true;
}

// Generic version for signed integers.
template <typename TInteger, typename TSource>
inline SEQAN_FUNC_ENABLE_IF(Is<SignedIntegerConcept<TInteger> >, bool)
lexicalCast(TInteger & target, TSource const & source)
{
    typedef typename Iterator<TSource const, Standard>::Type TIter;

    TIter it = begin(source, Standard());
    TIter itEnd = end(source, Standard());

    if (SEQAN_UNLIKELY(it == itEnd))
        return false;

    TInteger val = 0;

    if (*it != '-')
    {
        do
        {
            unsigned char digit = *it++ - '0';

            // invalid digit detection
            if (SEQAN_UNLIKELY(digit > 9))
                return false;

            // overflow detection
            if (SEQAN_UNLIKELY(val > MaxValue<TInteger>::VALUE / 10))
                return false;
            val *= 10;

            // overflow detection
            val += digit;
            if (SEQAN_UNLIKELY(val < digit))
                return false;
        }
        while (it != itEnd);
    }
    else
    {
        if (SEQAN_UNLIKELY(++it == itEnd))
            return false;
        do
        {
            unsigned char digit = *it++ - '0';

            // invalid digit detection
            if (SEQAN_UNLIKELY(digit > 9))
                return false;

            // overflow detection
            if (SEQAN_UNLIKELY(val < MinValue<TInteger>::VALUE / 10))
                return false;
            val *= 10;

            // overflow detection
            if (SEQAN_UNLIKELY(MinValue<TInteger>::VALUE - val > -(TInteger)digit))
                return false;
            val -= digit;
        }
        while (it != itEnd);
    }
    target = val;
    return true;
}

// Specialization for float.

template <typename TSource>
inline bool lexicalCast(float & target, TSource const & source)
{
    int offset;
    return (sscanf(toCString(source), "%g%n", &target, &offset) == 1) &&
           (static_cast<typename Size<TSource>::Type>(offset) == length(source));
}

// Specialization for double
template <typename TSource>
inline bool lexicalCast(double & target, TSource const & source)
{
    int offset;
    return (sscanf(toCString(source), "%lg%n", &target, &offset) == 1) &&
           (static_cast<typename Size<TSource>::Type>(offset) == length(source));
}

template <typename TTarget, typename TSource>
inline TTarget lexicalCast(TSource const & source)
{
    TTarget target;
    if (!lexicalCast(target, source))
        throw BadLexicalCast(target, source);
    return target;
}

// ----------------------------------------------------------------------------
// Function lexicalCastWithException()
// ----------------------------------------------------------------------------

/*!
 * @fn lexicalCastWithException
 * @headerfile <seqan/stream.h>
 * @brief Interpret the character sequence in <tt>source</tt> as the numeric value and write to <tt>target</tt>.
 *
 * @signature void lexicalCastWithException(target, source);
 *
 * @param[out] target  A numeric value to write to.
 * @param[in]  source  The source sequence of <tt>char</tt> to convert.
 *
 * @throw BadLexicalCast in the case that casting failed
 *
 * See @link lexicalCast @endlink for examples.
 */

template <typename TTarget, typename TSource>
inline void lexicalCastWithException(TTarget & target, TSource const & source)
{
    if (!lexicalCast(target, source))
        throw BadLexicalCast(target, source);
}

}

#endif //def SEQAN_STREAM_LEXICAL_CAST_H
