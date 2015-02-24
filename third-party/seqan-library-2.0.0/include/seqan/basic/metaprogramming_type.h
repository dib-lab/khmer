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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Metaprogramming for querying and modifying types.
//
// This header contains metafunctions for querying information about types
// and modifying it, such as querying for const-ness.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_TYPE_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_TYPE_H_

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

// ----------------------------------------------------------------------------
// Metafunction IsSameType
// ----------------------------------------------------------------------------

/*!
 * @mfn IsSameType
 * @headerfile <seqan/basic.h>
 * @brief Metaprogramming type comparison.
 *
 * @signature IsSameType<T1, T2>::Type;
 * @signature IsSameType<T1, T2>::VALUE;
 *
 * @tparam T1 Left-hand argument.
 * @tparam T2 Right-hand argument.
 *
 * @return Type  Either True or False, depending on whether T1 is the same type as T2.
 * @return VALUE The same as <tt>Type::VALUE</tt>.
 */

template <typename Type1, typename Type2>
struct IsSameType : False {};

template <typename Type1>
struct IsSameType<Type1, Type1> : True {};

// ----------------------------------------------------------------------------
// Metafunction MakeUnsigned
// ----------------------------------------------------------------------------

/*!
 * @mfn MakeUnsigned
 * @headerfile <seqan/basic.h>
 * @brief Convert an integral value into its unsigned variant.
 *
 * Returns T itself if T is not signed.
 *
 * @signature MakeUnsigned<T>::Type;
 *
 * @tparam T Input integral type.
 *
 * @return Type The unsigned version of T.
 *
 * @section Remarks
 *
 * The function is defined for all builtin integral types and with a fallback to that returns T if T is not a built-in
 * integral value.  You can specialize the metafunction for your custom types.
 */

template <typename T>
struct MakeUnsigned
{
    typedef
        typename If<typename IsSameType<T, __int8>::Type,       __uint8,
        typename If<typename IsSameType<T, char>::Type,         unsigned char,
        typename If<typename IsSameType<T, signed char>::Type,  unsigned char,
        typename If<typename IsSameType<T, signed short>::Type, unsigned short,
        typename If<typename IsSameType<T, signed int>::Type,   unsigned int,
        typename If<typename IsSameType<T, signed long>::Type,  unsigned long,
        typename If<typename IsSameType<T, __int64>::Type,      __uint64, T
        >::Type>::Type>::Type>::Type>::Type>::Type>::Type Type;
};

template <typename T>
struct MakeUnsigned<T const>
{
    typedef typename MakeUnsigned<T>::Type const Type;
};

// TODO(holtgrew): Internal metafunction unnecessary now?

template <typename T>
struct MakeUnsigned_ : MakeUnsigned<T> {};

// ----------------------------------------------------------------------------
// Metafunction MakeSigned
// ----------------------------------------------------------------------------

/*!
 * @mfn MakeSigned
 * @headerfile <seqan/basic.h>
 * @brief Convert an integral value into its signed variant.
 *
 * Returns T if T is already signed.
 *
 * @signature MakeSigned<T>::Type;
 *
 * @tparam T Input integral type.
 *
 * @return Type The signed version of T.
 *
 * @section Remarks
 *
 * The function is defined for all builtin integral types and with a fallback to that returns T if T is not a built-in
 * integral value.  You can specialize the metafunction for your custom types.
 */

template <typename T>
struct MakeSigned
{
    typedef
        typename If<typename IsSameType<T, char>::Type,           signed char,
        typename If<typename IsSameType<T, __int8>::Type,         __int8,
        typename If<typename IsSameType<T, unsigned char>::Type,  signed char,
        typename If<typename IsSameType<T, unsigned short>::Type, signed short,
        typename If<typename IsSameType<T, unsigned int>::Type,   signed int,
        typename If<typename IsSameType<T, unsigned long>::Type,  signed long,
        typename If<typename IsSameType<T, __uint64>::Type,       __int64, T
        >::Type>::Type>::Type>::Type>::Type>::Type>::Type Type;
};

template <typename T>
struct MakeSigned<T const>
{
    typedef typename MakeSigned<T>::Type const Type;
};

// TODO(holtgrew): Internal metafunction unnecessary now?

template <typename T>
struct MakeSigned_ : MakeSigned<T> {};

// ----------------------------------------------------------------------------
// Metafunction RemoveReference
// ----------------------------------------------------------------------------

/*!
 * @mfn RemoveReference
 * @headerfile <seqan/basic.h>
 * @brief Converts a (reference) type into the same type without a reference.
 *
 * @signature RemoveReference<T>::Type;
 *
 * @tparam T The input type.
 *
 * @return Type A corresponding non-reference type, e.g. <tt>int</tt> for <tt>T = &amp; int</tt>.
 */

#ifdef SEQAN_CXX11_STANDARD

template <typename T>
struct RemoveReference
{
    typedef typename std::remove_reference<T>::type Type;
};

#else

template <typename T>
struct RemoveReference
{
    typedef T Type;
};

template <typename T>
struct RemoveReference<T &> : RemoveReference<T> {};

#endif

// ----------------------------------------------------------------------------
// Metafunction RemoveReference
// ----------------------------------------------------------------------------

/*!
 * @mfn RemovePointer
 * @headerfile <seqan/basic.h>
 * @brief Converts a (pointer) type into the same type without a pointer.
 *
 * @signature RemovePointer<T>::Type;
 *
 * @tparam T The input type.
 *
 * @return Type A corresponding non-pointer type, e.g. <tt>int</tt> for <tt>T = *int</tt>.
 */

#ifdef SEQAN_CXX11_STANDARD

template <typename T>
struct RemovePointer
{
    typedef typename std::remove_pointer<T>::type Type;
};

#else

template <typename T>
struct RemovePointer
{
    typedef T Type;
};

template <typename T>
struct RemovePointer<T *>
{
    typedef T Type;
};

template <typename T>
struct RemovePointer<T * const>
{
    typedef T Type;
};

#endif

template <typename T>
struct IsPointer : False {};

template <typename T>
struct IsPointer<T *> : True {};

template <typename T>
struct IsPointer<T * const> : True {};

// ----------------------------------------------------------------------------
// Metafunction RemoveConst
// ----------------------------------------------------------------------------

/*!
 * @mfn RemoveConst
 * @headerfile <seqan/basic.h>
 * @brief Converts a (const) type into the corresponding non-const type.
 *
 * @signature RemoveConst<T>::Type;
 *
 * @tparam T Input type.
 *
 * @return Type A corresponding non-const type, e.g. <tt>int</tt> for <tt>T = const int</tt>.
 */

template <typename T>
struct RemoveConst
{
    typedef T Type;
};

template <typename T>
struct RemoveConst<T const> : public RemoveConst<T> {};

template <typename T>
struct RemoveConst<T &>
{
    typedef typename RemoveConst<T>::Type & Type;
};

/*
template <typename T>
struct RemoveConst<T const *>
{
    typedef typename RemoveConst<T>::Type * Type;
};
*/

template <typename T, size_t I>
struct RemoveConst<T const [I]>
{
    typedef T Type[I];
};

// TODO(holtgrew): Internal metafunction superflous?
template <typename T>
struct RemoveConst_ : RemoveConst<T> {};

// ----------------------------------------------------------------------------
// Metafunction CopyConst_
// ----------------------------------------------------------------------------

// TODO(holtgrew): Make public, document.

// returns TTo const, if TFrom is const, TTo otherwise

template <typename TFrom, typename TTo>
struct CopyConst_
{
    typedef TTo Type;
};

template <typename TFrom, typename TTo>
struct CopyConst_<TFrom const, TTo>
{
    typedef TTo const Type;
};

// ----------------------------------------------------------------------------
// Metafunction IsConst_
// ----------------------------------------------------------------------------

// TODO(holtgrew): Make public, complete documentation.

template <typename T>
struct IsConst_ : False
{};

template <typename T>
struct IsConst_<T const> : True
{};

// ----------------------------------------------------------------------------
// Metafunction ClassIdentifier_
// ----------------------------------------------------------------------------

// TODO(holtgrew): Make public, complete documentation or deletion candidate.

template <typename T>
struct ClassIdentifier_
{
    static inline void *
    getID()
    {
        static bool _id_dummy;
        return &_id_dummy;
    }
};

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_TYPE_H_
