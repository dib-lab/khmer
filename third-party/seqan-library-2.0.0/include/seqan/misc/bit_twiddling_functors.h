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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Implements basic functors.
// ==========================================================================

#ifndef INCLUDE_SEQAN_MISC_BIT_TWIDDLING_FUNCTORS_H_
#define INCLUDE_SEQAN_MISC_BIT_TWIDDLING_FUNCTORS_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Functor FunctorIdentity
// ----------------------------------------------------------------------------

struct FunctorIdentity
{
    template <typename TValue>
    inline TValue const & operator()(TValue const & val) const
    {
        return val;
    }
};

// ----------------------------------------------------------------------------
// Functor FunctorBitwiseAnd
// ----------------------------------------------------------------------------

struct FunctorBitwiseAnd
{
    template<typename TValue>
    inline TValue operator()(TValue const & valLhs, TValue const & valRhs) const
    {
        return valLhs & valRhs;
    }
};

// ----------------------------------------------------------------------------
// Functor FunctorBitwiseOr
// ----------------------------------------------------------------------------

struct FunctorBitwiseOr
{
    template <typename TValue>
    inline TValue operator()(TValue const & valLhs, TValue const & valRhs) const
    {
        return valLhs | valRhs;
    }
};

// ----------------------------------------------------------------------------
// Functor FunctorBitwiseXor
// ----------------------------------------------------------------------------

struct FunctorBitwiseXor
{
    template <typename TValue>
    inline TValue operator()(TValue const & valLhs, TValue const & valRhs) const
    {
        return valLhs ^ valRhs;
    }
};

// ----------------------------------------------------------------------------
// Functor FunctorBitwiseNot
// ----------------------------------------------------------------------------

struct FunctorBitwiseNot
{
    template <typename TValue>
    inline TValue operator()(TValue const & val) const
    {
        return ~val;
    }
};

// ----------------------------------------------------------------------------
// Functor FunctorNested
// ----------------------------------------------------------------------------

template <typename TBinaryFunctor, typename TUnaryFunctor1, typename TUnaryFunctor2>
struct FunctorNested
{
    template <typename TValue>
    inline TValue operator()(TValue const & lhs, TValue const & rhs) const
    {
        return TBinaryFunctor()(TUnaryFunctor1()(lhs), TUnaryFunctor2()(rhs));
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}

#endif // INCLUDE_SEQAN_MISC_BIT_TWIDDLING_FUNCTORS_H_
