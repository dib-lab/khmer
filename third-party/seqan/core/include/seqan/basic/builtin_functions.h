// ==========================================================================
//                            builtin_functions.h
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
// Documentation for built-in functions.
//
// This is used for documenting that certain global functions and operators
// are overridden for some classes.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_BUILTIN_FUNCTIONS_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_BUILTIN_FUNCTIONS_H_

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

/**
.Function.operator=
..cat:C++ built-in
..summary:C++ built-in Assignment operator.

.Function.operator+
..cat:C++ built-in
..summary:C++ built-in addition operator.

.Function.operator+ (unary)
..cat:C++ built-in
..summary:C++ built-in unary plus (integer promotion) operator.

.Function.operator-
..cat:C++ built-in
..summary:C++ built-in subtraction operator.

.Function.operator- (unary)
..cat:C++ built-in
..summary:C++ built-in unary minus (additive inverse) operator.

.Function.operator*
..cat:C++ built-in
..summary:C++ built-in multiplication operator.

.Function.operator/
..cat:C++ built-in
..summary:C++ built-in division operator.

.Function.operator%
..cat:C++ built-in
..summary:C++ built-in modulo operator.

.Function.operator++ (prefix)
..cat:C++ built-in
..summary:C++ built-in prefix increment operator.

.Function.operator++ (suffix)
..cat:C++ built-in
..summary:C++ built-in suffix increment operator.

.Function.operator-- (prefix)
..cat:C++ built-in
..summary:C++ built-in prefix decrement operator.

.Function.operator-- (suffix)
..cat:C++ built-in
..summary:C++ built-in suffix decrement operator.

.Function.operator==
..cat:C++ built-in
..summary:C++ built-in equal comparison operator.

.Function.operator!=
..cat:C++ built-in
..summary:C++ built-in inequal comparison operator.

.Function.operator>
..cat:C++ built-in
..summary:C++ built-in greater-than comparison operator.

.Function.operator<
..cat:C++ built-in
..summary:C++ built-in less-than comparison operator.

.Function.operator>=
..cat:C++ built-in
..summary:C++ built-in greather-than-or-equal comparison operator.

.Function.operator<=
..cat:C++ built-in
..summary:C++ built-in less-than-or-equal comparison operator.

.Function.operator!
..cat:C++ built-in
..summary:C++ built-in logical negation operator.

.Function.operator&&
..cat:C++ built-in
..summary:C++ built-in logical AND operator.

.DISABLED.Function.operator\pipe\pipe|operator||
..cat:C++ built-in
..summary:C++ built-in logical OR operator.

.Function.operator~
..cat:C++ built-in
..summary:C++ built-in bitwise NOT operator.

.Function.operator&
..cat:C++ built-in
..summary:C++ built-in bitwise AND operator.

.DISABLED.Function.operator\pipe|operator\pipe
..cat:C++ built-in
..summary:C++ built-in bitwise OR operator.

.Function.operator^
..cat:C++ built-in
..summary:C++ built-in bitwise XOR operator.

.Function.operator<<
..cat:C++ built-in
..summary:C++ built-in bitwise left shift operator.

.Function.operator<< (Stream)
..cat:C++ built-in
..summary:C++ built-in bitwise put-to/stream insertion operator.

.Function.operator>>
..cat:C++ built-in
..summary:C++ built-in bitwise right shift operator.

.Function.operator>> (Stream)
..cat:C++ built-in
..summary:C++ built-in bitwise get-from/stream extraction operator.

.Function.operator+=
..cat:C++ built-in
..summary:C++ built-in addition assignment operator.

.Function.operator-=
..cat:C++ built-in
..summary:C++ built-in subtraction assignment operator.

.Function.operator*=
..cat:C++ built-in
..summary:C++ built-in multiplication assignment operator.

.Function.operator/=
..cat:C++ built-in
..summary:C++ built-in division assignment operator.

.Function.operator%=
..cat:C++ built-in
..summary:C++ built-in modulo assignment operator.

.Function.operator&=
..cat:C++ built-in
..summary:C++ built-in bitwise AND assignment operator.

.DISABLED.Function.operator\pipe=|operator|=
..cat:C++ built-in
..summary:C++ built-in bitwise OR assignment operator.

.Function.operator^=
..cat:C++ built-in
..summary:C++ built-in bitwise XOR assignment operator.

.Function.operator<<=
..cat:C++ built-in
..summary:C++ built-in bitwise left shift assignment operator.

.Function.operator>>=
..cat:C++ built-in
..summary:C++ built-in bitwise right shift assignment operator.

.Function.operator[]
..cat:C++ built-in
..summary:C++ built-in array subscript operator.

.Function.operator* (indirection)
..cat:C++ built-in
..summary:C++ built-in indirection/object-pointed-to-by operator.

.Function.operator& (reference)
..cat:C++ built-in
..summary:C++ built-in reference/address-of operator.

.Function.operator->
..cat:C++ built-in
..summary:C++ built-in structure dereference operator.

.Function.operator->*
..cat:C++ built-in
..summary:C++ built-in member-pointed-to-by-b-of-object-pointed-to-by-a operator.

.Function.operator()
..cat:C++ built-in
..summary:C++ built-in function call operator.

.Function.operator,
..cat:C++ built-in
..summary:C++ built-in comma operator.

.Function.cast operator
..cat:C++ built-in
..summary:C++ built-in cast operator.

.Function.operator new
..cat:C++ built-in
..summary:C++ built-in allocation operator.

.Function.operator new[]
..cat:C++ built-in
..summary:C++ built-in array allocation operator operator.

.Function.operator delete
..cat:C++ built-in
..summary:C++ built-in deallocation operator.

.Function.operator delete[]
..cat:C++ built-in
..summary:C++ built-in array deallocation operator.
 */

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_BUILTIN_FUNCTIONS_H_
