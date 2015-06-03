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
// Align-specific metafunctions.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_ALIGN_METAFUNCTIONS_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_ALIGN_METAFUNCTIONS_H_

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
// Metafunction Cols
// ----------------------------------------------------------------------------

template <typename T>
struct Cols;

// ----------------------------------------------------------------------------
// Metafunction Col
// ----------------------------------------------------------------------------

/*!
 * @mfn Align#Col
 * @headerfile <seqan/align.h>
 * @brief The column type for @link Align @endlink objects.
 *
 * @signature Col<TAlign>::Type
 *
 * @tparam TAlign The @link Align @endlink object to query for its column type.
 *
 * @return Type   The resulting type.
 */

template <typename T>
struct Col : Value<typename Cols<T>::Type>
{};

// ----------------------------------------------------------------------------
// Metafunction Rows
// ----------------------------------------------------------------------------

template <typename T>
struct Rows;

// ----------------------------------------------------------------------------
// Metafunction Row
// ----------------------------------------------------------------------------

template <typename T>
struct Row : Value<typename Rows<T>::Type>
{};

template <typename T>
struct Row<T const>
{
    typedef typename Row<T>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction StringSetType
// ----------------------------------------------------------------------------

template <typename T>
struct StringSetType;

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_ALIGN_METAFUNCTIONS_H_
