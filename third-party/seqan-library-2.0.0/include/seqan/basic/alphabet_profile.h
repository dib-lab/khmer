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
// Author: Tobias Rausch <rausch@embl.de>
// ==========================================================================
// Profile alphabet character code.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_ALPHABET_PROFILE_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_ALPHABET_PROFILE_H_

#include <seqan/misc/memset.h>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TSpec>
class Proxy;

template <typename TSpec>
struct IteratorProxy;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class ProfileChar
 *
 * @headerfile <seqan/basic.h>
 *
 * @brief Alphabet type for profiles over another alphabet.
 *
 * @signature template <typename TValue[, typename TCount[, typename TSpec]]>
 *            class ProfileChar;
 *
 * @tparam TValue The underlying alphabet type.
 * @tparam TCount The type to use for counting, default: <tt>unsigned int</tt>.
 * @tparam TSpec  Specialization tag, default: <tt>void</tt>
 */

/*!
 * @var VariableType ProfileChar::count[]
 *
 * @brief Array of ValueSize elements, giving counts in profile.
 */

template <typename TValue, typename TCount = unsigned, typename TSpec = void>
class ProfileChar;

template <typename TValue, typename TCount, typename TSpec>
class ProfileChar
{
public:
    typedef typename ValueSize<ProfileChar>::Type TSize;

    TCount count[ValueSize<ProfileChar>::VALUE];

    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    ProfileChar()
    {
        memset<ValueSize<ProfileChar>::VALUE * sizeof(TCount), (unsigned char) 0>(count);
    }

    ProfileChar(ProfileChar const & other_data)
    {
        for (TSize i = 0; i < ValueSize<ProfileChar>::VALUE; ++i)
            count[i] = other_data.count[i];
    }

    // TODO(holtgrew): Limit TOther to SourceValue?
    template <typename TOther>
    ProfileChar(TOther const & other_data)
    {
        memset<ValueSize<ProfileChar>::VALUE * sizeof(TCount), (unsigned char) 0>(count);
        count[ordValue(TValue(other_data))] = 1;
    }

    template <typename TSpec2>
    ProfileChar(Proxy<TSpec2> const & proxy)
    {
        memset<ValueSize<ProfileChar>::VALUE * sizeof(TCount), (unsigned char) 0>(count);
        assign(*this, getValue(proxy));
    }

    // ------------------------------------------------------------------------
    // Assignment operators;  Have to be defined in class.
    // ------------------------------------------------------------------------

    ProfileChar &
    operator=(ProfileChar const & other_data)
    {
        if (this == &other_data) return *this;

        for (TSize i = 0; i < ValueSize<ProfileChar>::VALUE; ++i)
            count[i] = other_data.count[i];
        return *this;
    }

    template <typename TOther>
    ProfileChar &
    operator=(TOther const & other_data)
    {
        memset<ValueSize<ProfileChar>::VALUE * sizeof(TCount), 0u>(count);
        count[ordValue(TValue(other_data))] = 1;
        return *this;
    }

    // ------------------------------------------------------------------------
    // Type conversion operators;  Have to be defined in class.
    // ------------------------------------------------------------------------

    operator char()
    {
        typename Size<ProfileChar>::Type maxIndex = getMaxIndex(*this);
        return (maxIndex == ValueSize<ProfileChar>::VALUE - 1) ? gapValue<char>() : (char) TValue(maxIndex);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ValueSize
// ----------------------------------------------------------------------------

/*!
 * @mfn ProfileChar#ValueSize
 * @brief Number of different values a value type object can have.
 *
 * @signature ValueSize<T>::VALUE;
 *
 * @tparam T The type to query.
 *
 * @return VALUE Number of different values T can have.
 */

template <typename TValue, typename TCount, typename TSpec>
struct ValueSize<ProfileChar<TValue, TCount, TSpec> >
{
    enum { VALUE = ValueSize<TValue>::VALUE + 1 };
    typedef unsigned Type;
};

// ----------------------------------------------------------------------------
// Metafunction SourceValue
// ----------------------------------------------------------------------------

/*!
 * @mfn ProfileChar#SourceValue
 * @brief Returns underlying value for ProfileChar.
 *
 * @signature SourceValue<T>::Type
 *
 * @tparam T Type to query.
 *
 * @return Type The type of the underlying character.
 *
 * @section Examples
 *
 * @code{.cpp}
 * typedef ProfileChar<Dna5>               TProfileChar;
 * typedef SourceValue<TProfileChar>::Type TType;  // Is Dna.
 * @endcode
 */

template <typename T>
struct SourceValue;

template <typename TValue, typename TCount, typename TSpec>
struct SourceValue<ProfileChar<TValue, TCount, TSpec> >
{
    typedef TValue Type;
};

template <typename TValue, typename TCount, typename TSpec>
struct SourceValue<ProfileChar<TValue, TCount, TSpec> const> :
            SourceValue<ProfileChar<TValue, TCount, TSpec> >
{};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function operator==()
// ----------------------------------------------------------------------------

template <typename TValue, typename TCount, typename TSpec>
inline bool
operator==(ProfileChar<TValue, TCount, TSpec> const & lhs,
           ProfileChar<TValue, TCount, TSpec> const & rhs)
{
    typedef ProfileChar<TValue, TCount, TSpec> TProfileChar;
    typedef typename ValueSize<TProfileChar>::Type TValueSize;

    for (TValueSize i = 0; i < ValueSize<TProfileChar>::VALUE; ++i)
        if (lhs.count[i] != rhs.count[i])
            return false;
    return true;
}

// ----------------------------------------------------------------------------
// Function operator!=()
// ----------------------------------------------------------------------------

template <typename TValue, typename TCount, typename TSpec>
inline
bool
operator!=(ProfileChar<TValue, TCount, TSpec> const & lhs,
           ProfileChar<TValue, TCount, TSpec> const & rhs)
{
    typedef ProfileChar<TValue, TCount, TSpec> TProfileChar;
    typedef typename ValueSize<TProfileChar>::Type TSize;

    for (TSize i = 0; i < ValueSize<TProfileChar>::VALUE; ++i)
        if (lhs.count[i] != rhs.count[i])
            return true;
    return false;
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

/*!
 * @fn ProfileChar#empty
 * @brief Check whether there are only gaps in the representation of the ProfileChar.
 *
 * @signature bool empty(c);
 *
 * @param  c    ProfileChar to query.
 * @return bool Whether or not the ProfileChar only contains gaps.
 */
template <typename TSourceValue, typename TSourceCount, typename TSourceSpec>
bool empty(ProfileChar<TSourceValue, TSourceCount, TSourceSpec> const & source)
{
    typedef typename ValueSize<ProfileChar<TSourceValue, TSourceCount, TSourceSpec> const>::Type TSize;

    for (TSize i = 0; i < ValueSize<TSourceValue>::VALUE; ++i)
        if (source.count[i])
            return false;
    return true;
}

// ----------------------------------------------------------------------------
// Helper Function getMaxIndex()
// ----------------------------------------------------------------------------

/*!
 * @fn ProfileChar#getMaxIndex
 * @brief Return number of dominating entry in ProfileChar.
 *
 * @signature TSize getMaxIndex(c);
 *
 * @param[in] c     ProfileChar to query for its dominating entry.
 * @return    TSize index (with the @link FiniteOrderedAlphabetConcept#ordValue @endlink) of the dominating character
 *                  in <tt>c</tt>
 */
template <typename TSourceValue, typename TSourceCount, typename TSourceSpec>
typename Size<ProfileChar<TSourceValue, TSourceCount, TSourceSpec> const>::Type
getMaxIndex(ProfileChar<TSourceValue, TSourceCount, TSourceSpec> const & source)
{
    typedef ProfileChar<TSourceValue, TSourceCount, TSourceSpec> TProfileChar;
    typedef typename Size<TProfileChar>::Type TSize;
    TSize maxIndex = 0;
    TSourceCount maxCount = source.count[0];
    for (TSize i = 1; i < ValueSize<TProfileChar>::VALUE; ++i)
    {
        if (source.count[i] > maxCount)
        {
            maxIndex = i;
            maxCount = source.count[i];
        }
    }
    return maxIndex;
}

// ----------------------------------------------------------------------------
// Helper Function totalCount()
// ----------------------------------------------------------------------------

/*!
 * @fn ProfileChar#totalCount
 * @brief Return sum of counts in ProfileChar.
 *
 * @signature TCount totalCount(c);
 *
 * @param[in] c      ProfileChar to query.
 * @return    TCount Total number of characters represented by <tt>c</tt>.
 */
template <typename TSourceValue, typename TSourceCount, typename TSourceSpec>
TSourceCount totalCount(ProfileChar<TSourceValue, TSourceCount, TSourceSpec> const & source)
{
    typedef ProfileChar<TSourceValue, TSourceCount, TSourceSpec> TProfileChar;
    typedef typename Size<TProfileChar>::Type TSize;
    TSourceCount totalCount = source.count[0];
    for (TSize i = 1; i < ValueSize<TProfileChar>::VALUE; ++i)
        totalCount += source.count[i];
    return totalCount;
}


// ----------------------------------------------------------------------------
// Function assign()
// ----------------------------------------------------------------------------

// TODO(holtgrew): What if there only are gaps?

template <typename TTargetValue, typename TTargetSpec, typename TSourceValue, typename TSourceCount, typename TSourceSpec>
inline void
assign(SimpleType<TTargetValue, TTargetSpec> & target,
       ProfileChar<TSourceValue, TSourceCount, TSourceSpec> const & source)
{
    target.value = getMaxIndex(source);
}

// ----------------------------------------------------------------------------
// Function convertImpl()
// ----------------------------------------------------------------------------

template <typename TTarget, typename T, typename TSourceValue, typename TSourceCount, typename TSourceSpec>
inline typename Convert<TTarget, ProfileChar<TSourceValue, TSourceCount, TSourceSpec> >::Type
convertImpl(Convert<TTarget, T> const &,
            ProfileChar<TSourceValue, TSourceCount, TSourceSpec> const & source)
{
    return (getMaxIndex(source) == ValueSize<TSourceValue>::VALUE) ? convertImpl(Convert<TTarget, T>(), '-') : convertImpl(Convert<TTarget, T>(), TSourceValue(getMaxIndex(source)));
}

// ----------------------------------------------------------------------------
// Function operator<<();  Stream output.
// ----------------------------------------------------------------------------

template <typename TStream, typename TValue, typename TCount, typename TSpec>
inline TStream &
operator<<(TStream & os, ProfileChar<TValue, TCount, TSpec> const & rhs)
{
    typedef ProfileChar<TValue, TCount, TSpec> TProfileChar;
    typedef typename Size<TProfileChar>::Type TSize;
    for (TSize i = 0; i < ValueSize<TProfileChar>::VALUE; ++i)
        os << i << ':' << rhs.count[i] << ' ' << ';';
    return os;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_ALPHABET_PROFILE_H_
