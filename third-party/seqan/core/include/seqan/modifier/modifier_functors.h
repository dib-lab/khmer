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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_MODIFIER_MODIFIER_FUNCTORS_H_
#define SEQAN_MODIFIER_MODIFIER_FUNCTORS_H_

#include <cctype>

// TODO(holtgrew): Make the structs here into classes.

namespace seqan
{

// ==========================================================================
// Forwards
// ==========================================================================

// ==========================================================================
// Classes, Enums, Typedefs
// ==========================================================================

// --------------------------------------------------------------------------
// Class FunctorUpcase
// --------------------------------------------------------------------------

/**
.Class.FunctorUpcase
..cat:Modifier
..summary:Functor that returns the upper case character to a given character.
..signature:FunctorUpcase<TValue>
..param.TValue:The input value type.
..remarks:This Functor is a derivation of the STL unary function.
..include:seqan/modifier.h
*/

template <typename InType, typename Result = InType>
struct FunctorUpcase : public std::unary_function<InType, Result> 
{
    inline Result operator()(InType x) const
    {
        return toupper(x);
    }
};

// --------------------------------------------------------------------------
// Class FunctorLowcase
// --------------------------------------------------------------------------

/**
.Class.FunctorLowcase
..cat:Modifier
..summary:Functor that returns the lower case character to a given character.
..signature:FunctorLowcase<TValue>
..param.TValue:The input value type.
..remarks:This Functor is a derivation of the STL unary function.
..include:seqan/modifier.h
*/

template <typename InType, typename Result = InType>
struct FunctorLowcase : public std::unary_function<InType, Result> 
{
    inline Result operator()(InType x) const
    {
        return tolower(x);
    }
};

// --------------------------------------------------------------------------
// Class FunctorConvert
// --------------------------------------------------------------------------

/**
.Class.FunctorConvert
..cat:Modifier
..summary:Functor that converts a $TInValue$ type to a $TOutValue$ type character.
..signature:FunctorConvert<TInValue, TOutValue>
..param.TInValue:The input value type.
..param.TOutValue:The output value type.
..remarks:This Functor is a derivation of the STL unary function.
..include:seqan/modifier.h
*/

template <typename InType, typename OutType>
struct FunctorConvert : public std::unary_function<InType,OutType> 
{
    inline OutType operator()(InType x) const
    {
        return x; 
    }
};

// --------------------------------------------------------------------------
// Helper Structs with Translation Tables
// --------------------------------------------------------------------------

// Manual forward for the complementing of characters.
template <typename TValue> struct FunctorComplement;


template <typename T = void>
struct TranslateTableDna5ToDna5Complement_
{
    static char const VALUE[5];
};

template <typename T>
char const TranslateTableDna5ToDna5Complement_<T>::VALUE[5] = {'T', 'G', 'C', 'A', 'N'};

template <typename T = void>
struct TranslateTableRna5ToRna5Complement_
{
    static char const VALUE[5];
};


template <typename T>
char const TranslateTableRna5ToRna5Complement_<T>::VALUE[5] = {'U', 'G', 'C', 'A', 'N'};

// --------------------------------------------------------------------------
// Class FunctorComplement
// --------------------------------------------------------------------------

/**
.Class.FunctorComplement:
..cat:Modifier
..summary:Functor that returns the complement nucleotide to a given nucleotide.
..signature:FunctorComplement<TValue>
..param.TValue:The input value type.
...type:Spec.Dna
...type:Spec.Dna5
...type:Spec.Rna
...type:Spec.Rna5
..remarks:This Functor is a derivation of the STL unary function.
..include:seqan/modifier.h
*/

template <>
struct FunctorComplement<char> : public std::unary_function<Dna5,Dna5> 
{
    inline Dna5 operator()(Dna5 x) const
    {
        return TranslateTableDna5ToDna5Complement_<>::VALUE[x.value]; 
    }
};

template <>
struct FunctorComplement<Dna> : public std::unary_function<Dna,Dna> 
{
    inline Dna operator()(Dna x) const
    {
        return TranslateTableDna5ToDna5Complement_<>::VALUE[x.value]; 
    }
};


template <>
struct FunctorComplement<Dna5> : public std::unary_function<Dna5,Dna5> 
{
    inline Dna5 operator()(Dna5 x) const
    {
        return TranslateTableDna5ToDna5Complement_<>::VALUE[x.value]; 
    }
};


template <>
struct FunctorComplement<Rna> : public std::unary_function<Rna,Rna> 
{
    inline Rna operator()(Rna x) const
    {
        return TranslateTableRna5ToRna5Complement_<>::VALUE[x.value]; 
    }
};


template <>
struct FunctorComplement<Rna5> : public std::unary_function<Rna5,Rna5> 
{
    inline Dna5 operator()(Rna5 x) const
    {
        return TranslateTableRna5ToRna5Complement_<>::VALUE[x.value]; 
    }
};

template <>
struct FunctorComplement<DnaQ> : public std::unary_function<DnaQ,DnaQ> 
{
    inline DnaQ operator()(DnaQ x) const
    {
        int qual = getQualityValue(x);
        x = TranslateTableDna5ToDna5Complement_<>::VALUE[ordValue((Dna)x)];
        assignQualityValue(x, qual);
        return x;
    }
};

template <>
struct FunctorComplement<Dna5Q> : public std::unary_function<Dna5Q,Dna5Q> 
{
    inline Dna5Q operator()(Dna5Q x) const {
        int qual = getQualityValue(x);
        x = TranslateTableDna5ToDna5Complement_<>::VALUE[ordValue((Dna5)x)];
        assignQualityValue(x, qual);
        return x;
    }
};

}  // namespace seqan

#endif  // SEQAN_MODIFIER_MODIFIER_FUNCTORS_H_
