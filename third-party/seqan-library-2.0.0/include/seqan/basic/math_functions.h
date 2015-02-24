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
// Author: Andreas Gogol-Döring <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Math-related utility routines.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_MATH_FUNCTIONS_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_MATH_FUNCTIONS_H_

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

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _intPow()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Document and make public.

template <typename TValue, typename TExponent>
inline TValue _intPow(TValue a, TExponent b)
{
    SEQAN_CHECKPOINT;
    TValue ret = 1;
    while (b != 0) {
        if (b & 1) ret *= a;
        a *= a;
        b >>= 1;
    }
    return ret;
}

// ----------------------------------------------------------------------------
// Function log2()
// ----------------------------------------------------------------------------

/*!
 * @fn log2
 * @headerfile <seqan/basic.h>
 * @brief Computes floored logarithm of base 2 for integer types
 *
 * @signature unsigned log2(i);
 *
 * @param[in] i An integer type.
 *
 * @return unsigned The largest integer smaller or equal than the logarithm of <tt>i</tt>.
 */

// TODO(holtgrew): Should this maybe called log2floor for consistency with Log2Floor<>::VALUE?

template <int BITS_MAX>
struct Log2Impl_
{
    template <typename T>
    static inline unsigned int
    log2(T val, unsigned int offset)
    {
        unsigned int val2 = val >> (BITS_MAX / 2);
        if (val2)
        {
            val = val2;
            offset += BITS_MAX / 2;
        }
        return Log2Impl_<BITS_MAX / 2>::log2(val, offset);
    }
};

template <>
struct Log2Impl_<1>
{
    template <typename T>
    static inline unsigned int
    log2(T /*val*/, unsigned int offset)
    {
        return offset;
    }
};

template <typename T>
inline unsigned int
log2(T val)
{
    enum
    {
//      BITS_PER_VALUE = BitsPerValue<T>::VALUE //  TODO(holtgrew): portable bits-per-char!
        BITS_PER_VALUE = sizeof(T) * 8
    };

    return Log2Impl_<BITS_PER_VALUE>::log2(val, 0);
}

// ----------------------------------------------------------------------------
// Function _min()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Subject to removal.  http://trac.mi.fu-berlin.de/seqan/ticket/855

// to avoid conflicts with non-standard macros and namespaces
// we define our own Min/Max functions

template<typename Tx_>
SEQAN_HOST_DEVICE inline
const Tx_& _min(const Tx_& _Left, const Tx_& Right_)
{   // return smaller of _Left and Right_
    if (_Left < Right_)
        return _Left;
    else
        return Right_;
}

template<typename Tx_, typename Ty_>
SEQAN_HOST_DEVICE inline
Tx_ _min(const Tx_& _Left, const Ty_& Right_)
{   // return smaller of _Left and Right_
    return (Right_ < _Left ? Right_ : _Left);
}

// ----------------------------------------------------------------------------
// Function _max()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Subject to removal.  http://trac.mi.fu-berlin.de/seqan/ticket/855

// to avoid conflicts with non-standard macros and namespaces
// we define our own Min/Max functions

template<typename Ty_>
inline Ty_ const &
_max(const Ty_& _Left, const Ty_& Right_)
{   // return larger of _Left and Right_
    if (_Left < Right_)
        return Right_;
    else
        return _Left;
}

template<typename Tx_, typename Ty_>
inline Tx_
_max(const Tx_& _Left, const Ty_& Right_)
{   // return smaller of _Left and Right_
    return (Right_ < _Left ? _Left : Right_);
}

// ----------------------------------------------------------------------------
// Function _abs()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Make public, document.  This is here since cmath's abs is only defined for floats/doubles.

template <typename T>
inline
T _abs(T const & x)
{
    if (x < static_cast<T>(0))
        return -x;
    else
        return x;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_MATH_FUNCTIONS_H_

