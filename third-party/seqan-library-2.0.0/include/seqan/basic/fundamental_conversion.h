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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Fundamental conversion code.
// ==========================================================================

// TODO(holtgrew): Move to alphabet submodule?

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_FUNDAMENTAL_CONVERSION_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_FUNDAMENTAL_CONVERSION_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

/*!
 * @mfn Convert
 * @headerfile <seqan/basic.h>
 * @brief Return type of a conversion.
 *
 * @signature Convert<TTarget, TSource>::Type;
 *
 * @tparam TSource Type of the object that should be converted to <tt>Target</tt>.
 * @tparam TTarget Type the object should be converted to.
 *
 * If instances of <tt>Source</tt> can be re-interpreted as instances of <tt>Target</tt>, than this metafunction returns
 * a reference, otherwise it returns <tt>Target</tt>, that is convert returns a temporary.
 *
 * @section Remarks
 *
 * A constant instance of <tt>Convert</tt> is (ab)used as tag argument of convertImpl.
 */

template <typename TTarget, typename TSource = void>
struct Convert
{
    typedef TTarget Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function convertImpl()
// ----------------------------------------------------------------------------

/*!
 * @fn convertImpl
 * @headerfile <seqan/basic.h>
 * @brief Implements convert.
 *
 * @signature T convertImpl(convert, source);
 *
 * @param[in] convert Object that specifies the conversion.  A constant instance of Convert is used to specify the
 *                    conversion target.
 * @param[in] source  An object that should be converted.
 *
 * @return T <tt>source</tt> converted to the type specified by convert.
 *
 * @section Remarks
 *
 * This function implements convert.  It is recommended to use convert rather than <tt>convertImpl</tt>.
 *
 * @see convert
 */

// NOTE(doering): Specialize convertImpl, use convert.
// NOTE(doering): Conversion of one char into another char happes at another place.
// NOTE(doering): Conversion of sequences happens at another place.
// NOTE(doering): Can copy or reinterpret, depending on Convert::Type

template <typename TTarget, typename T, typename TSource>
SEQAN_HOST_DEVICE inline typename Convert<TTarget, TSource>::Type
convertImpl(Convert<TTarget, T> const,
            TSource const & source)
{
    return source;
}

// ----------------------------------------------------------------------------
// Function convert()
// ----------------------------------------------------------------------------

/*!
 * @fn convert
 * @headerfile <seqan/basic.h>
 * @brief Converts a value into another value.
 *
 * @signature T convert<Target>(source);
 *
 * @tparam Target The type <tt>source</tt> is converted to.
 *
 * @param[in] source An object that is converted to <tt>Target</tt>.
 *
 * @return TReturn <tt>source</tt> converted to <tt>Target</tt>.  If <tt>source</tt> can be re-interpreted as instance
 *                 of <tt>Target</tt>, then a reference is returned.  Otherwise the function returns a temporary
 *                 object.
 *
 * This function is implemented in convertImpl. Do not specialize <tt>convert</tt>, specialize convertImpl instead.
 *
 * @see convertImpl
 */

template <typename TTarget, typename TSource>
SEQAN_HOST_DEVICE inline typename Convert<TTarget, TSource>::Type
convert(TSource const & source)
{
    return convertImpl(Convert<TTarget, TSource>(), source);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_FUNDAMENTAL_CONVERSION_H_
