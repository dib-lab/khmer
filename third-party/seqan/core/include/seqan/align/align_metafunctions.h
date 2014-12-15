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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Align-specific metafunctions.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_ALIGN_METAFUNCTIONS_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_ALIGN_METAFUNCTIONS_H_

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

/**
.Metafunction.Cols:
..cat:Alignments
..summary:Type of column container of an alignment.
..signature:Cols<T>::Type
..param.T:An alignment type.
...type:Class.Align
..returns.param.Type:The type of the container that allows access to the columns of $T$.
..include:seqan/align.h
*/

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

/**
.Metafunction.Col:
..cat:Alignments
..summary:Type of a column in an alignment.
..signature:Col<T>::Type
..param.T:An alignment type.
...type:Class.Align
..returns.param.Type:The column type of $T$.
..remarks:The returned type is equivalent to $Value<Cols<T>::Type>::Type$.
..see:Metafunction.Cols
..see:Metafunction.Value
..include:seqan/align.h
*/

template <typename T>
struct Col : Value<typename Cols<T>::Type>
{};

// ----------------------------------------------------------------------------
// Metafunction Rows
// ----------------------------------------------------------------------------

/**
.Metafunction.Rows:
..cat:Alignments
..summary:Type of row container of an alignment.
..signature:Rows<T>::Type
..param.T:An alignment type.
...type:Class.Align
..returns.param.Type:The type of the container that allows access to the rows of $T$.
..see:Metafunction.Cols
..include:seqan/align.h
*/

template <typename T>
struct Rows;

// ----------------------------------------------------------------------------
// Metafunction Row
// ----------------------------------------------------------------------------

/**
.Metafunction.Row:
..cat:Alignments
..summary:Type of a row in an alignment.
..signature:Row<T>::Type
..param.T:An alignment type.
...type:Class.Align
..returns.param.Type:The row type of $T$.
..remarks:The returned type is equivalent to $Value<Rows<T>::Type>::Type$.
..see:Metafunction.Rows
..see:Metafunction.Value
..include:seqan/align.h
*/

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

/**
.Metafunction.StringSetType:
..cat:Alignments
..summary:Return type of @Function.stringSet@ function.
..signature:StringSetType<T>::Type
..param.T:Alignment data structure.
..param.T.type:Spec.Alignment Graph
..param.T.type:Class.Align
..returns.param.Type:A @Class.StringSet.string set@ type of a reference to a string set type.
..see:Function.stringSet
..include:seqan/align.h
*/
template <typename T>
struct StringSetType;

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_ALIGN_METAFUNCTIONS_H_
