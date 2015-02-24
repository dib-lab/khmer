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
// Author: Hannes Hauswedell <hannes.hauswedell@fu-berlin.de>
// ==========================================================================
// Reduced Versions of the 24-letter amino acid alphabet
// ==========================================================================

#ifndef SEQAN_REDUCED_AMINOACID_BASE_H_
#define SEQAN_REDUCED_AMINOACID_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// -----------------------------------------------------------------------
// Tag ReducedAminoAcid
// -----------------------------------------------------------------------

/*!
 * @class ReducedAminoAcid
 * @extends SimpleType
 * @brief Reduced versions of the amino acid alphabet.
 * @headerfile seqan/reduced_aminoacid.h
 *
 * @signature template <typename TReductionSpec>
 * using ReducedAminoAcid = SimpleType<unsigned char, ReducedAminoAcid_<TReductionSpec> >;
 *
 * @tparam TReductionSpec Currently only @link Murphy10 @endlink available
 *
 * @section Remarks
 *
 * The alias template is only available when SEQAN_CXX11_STANDARD is defined
 * and your compiler supports alias templates (Visual Studio >= 2006-2014, any fairly
 * recent Clang, GCC). Otherwise you have to use the underscored type and
 * the full definition, i.e.
 * <tt>SimpleType&lt;unsigned char, ReducedAminoAcid_&lt;TReductionSpec&gt; &gt;</tt>.
 *
 * @see Murphy10
 */

template <typename TRedSpec>
struct ReducedAminoAcid_ {};

#if defined (SEQAN_CXX11_STANDARD) && ( !defined (_MSC_VER) || _MSC_VER >= 1800 )
template <typename TRedSpec>
using ReducedAminoAcid = SimpleType<unsigned char, ReducedAminoAcid_<TRedSpec> >;
#endif

// ============================================================================
// Metafunctions
// ============================================================================

// -----------------------------------------------------------------------
// Metafunction CompareType
// -----------------------------------------------------------------------

template <typename TRedSpec>
struct CompareType<SimpleType<unsigned char, ReducedAminoAcid_<TRedSpec> >,
                   __uint8>
{
    typedef SimpleType<unsigned char, ReducedAminoAcid_<TRedSpec> > Type;
};

template <typename TRedSpec>
struct CompareType<SimpleType<unsigned char, ReducedAminoAcid_<TRedSpec> >,
                   char>
{
    typedef SimpleType<unsigned char, ReducedAminoAcid_<TRedSpec> > Type;
};

template <typename TRedSpec>
struct CompareType<SimpleType<unsigned char, ReducedAminoAcid_<TRedSpec> >,
                   AminoAcid>
{
    typedef SimpleType<unsigned char, ReducedAminoAcid_<TRedSpec> > Type;
};

// -----------------------------------------------------------------------
// Translation Tables (implementations see extra files)
// -----------------------------------------------------------------------

template <typename TRedSpec, typename TSpec = void>
struct TranslateTableCharToRedAA_;

template <typename TRedSpec, typename TSpec = void>
struct TranslateTableAAToRedAA_;

template <typename TRedSpec, typename TSpec = void>
struct TranslateTableByteToRedAA_;

template <typename TRedSpec, typename TSpec = void>
struct TranslateTableRedAAToChar_;

// ============================================================================
// Functions
// ============================================================================

}
#endif // def SEQAN_REDUCED_AMINOACID_BASE_H_
