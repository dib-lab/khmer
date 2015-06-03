// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
// Copyright (c) 2013 NVIDIA Corporation
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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================

//SEQAN_NO_DDDOC:do not generate documentation for this file

#ifndef INDEX_FM_COMPRESSED_SA_ITERATOR_H_
#define INDEX_FM_COMPRESSED_SA_ITERATOR_H_

namespace seqan {

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig>
struct Iterator<CompressedSA<TText, TSpec, TConfig> const, Standard>
{
    typedef Iter<CompressedSA<TText, TSpec, TConfig> const, PositionIterator> Type;
};

template <typename TText, typename TSpec, typename TConfig>
struct Iterator<CompressedSA<TText, TSpec, TConfig>, Standard>
{
    typedef Iter<CompressedSA<TText, TSpec, TConfig>, PositionIterator> Type;
};

template <typename TText, typename TSpec, typename TConfig>
struct Iterator<CompressedSA<TText, TSpec, TConfig>, Rooted>:
    Iterator<CompressedSA<TText, TSpec, TConfig>, Standard>{};

template <typename TText, typename TSpec, typename TConfig>
struct Iterator<CompressedSA<TText, TSpec, TConfig> const, Rooted>:
    Iterator<CompressedSA<TText, TSpec, TConfig> const, Standard>{};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------
/*!
 * @fn CompressedSA#begin
 * @headerfile <seqan/index.h>
 * @brief Returns an iterator pointing to the first position of a compresses suffix array.
 *
 * @signature TIterator begin(compressedSA, tag);
 *
 * @param[in] compressedSA The compresses suffix array to be traversed.
 * @param[in] tag The specialisation of the iterator to be returned by the function. Types: @link
 *                ContainerIteratorTags#Standard @endlink, @link ContainerIteratorTags#Rooted @endlink
 *
 * @return TIterator Returns an iterator pointing to the first position of a compresses suffix array.  Types: <tt>The
 *                   result of Iterator&lt;Index&lt;TText, TIndexSpec&gt;, TSpec&gt;::Type</tt>
 */
template <typename TText, typename TSpec, typename TConfig>
inline typename Iterator<CompressedSA<TText, TSpec, TConfig>, Standard>::Type
begin(CompressedSA<TText, TSpec, TConfig> & compressedSA, Standard const & /* dummy */)
{
    return typename Iterator<CompressedSA<TText, TSpec, TConfig>, Standard>::Type(compressedSA, 0);
}

template <typename TText, typename TSpec, typename TConfig>
inline typename Iterator<CompressedSA<TText, TSpec, TConfig> const, Standard>::Type
begin(CompressedSA<TText, TSpec, TConfig> const & compressedSA, Standard const & /* dummy */)
{
    return typename Iterator<CompressedSA<TText, TSpec, TConfig> const, Standard>::Type(compressedSA, 0);
}

template <typename TText, typename TSpec, typename TConfig>
inline typename Iterator<CompressedSA<TText, TSpec, TConfig>, Rooted>::Type
begin(CompressedSA<TText, TSpec, TConfig> & compressedSA, Rooted const & /* dummy */)
{
    return typename Iterator<CompressedSA<TText, TSpec, TConfig>, Rooted>::Type(compressedSA, 0);
}

template <typename TText, typename TSpec, typename TConfig>
inline typename Iterator<CompressedSA<TText, TSpec, TConfig> const, Rooted>::Type
begin(CompressedSA<TText, TSpec, TConfig> const & compressedSA, Rooted const & /* dummy */)
{
    return typename Iterator<CompressedSA<TText, TSpec, TConfig> const, Rooted>::Type(compressedSA, 0);
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------
/*!
 * @fn CompressedSA#end
 * @headerfile <seqan/index.h>
 * @brief Returns an iterator pointing to the position behind the last element of a compresses suffix array.
 *
 * @signature TIterator end(compressedSA, tag);
 *
 * @param[in] compressedSA The compresses suffix array to be traversed.
 * @param[in] tag The specialisation of the iterator to be returned by the function. Types: @link
 *                ContainerIteratorTags#Standard @endlink, @link ContainerIteratorTags#Rooted @endlink
 *
 * @return TIterator Returns an iterator pointing to the position behind the last element of a compresses suffix array.
 *                   Types: <tt>The result of Iterator&lt;Index&lt;TText, TIndexSpec&gt;, TSpec&gt;::Type</tt>
 */

template <typename TText, typename TSpec, typename TConfig>
inline typename Iterator<CompressedSA<TText, TSpec, TConfig>, Rooted>::Type
end(CompressedSA<TText, TSpec, TConfig> & compressedSA, Rooted const & /* dummy */)
{
    return typename Iterator<CompressedSA<TText, TSpec, TConfig>, Rooted>::Type(compressedSA, length(compressedSA));
}

template <typename TText, typename TSpec, typename TConfig>
inline typename Iterator<CompressedSA<TText, TSpec, TConfig> const, Rooted>::Type
end(CompressedSA<TText, TSpec, TConfig> const & compressedSA, Rooted/* dummy */)
{
    return typename Iterator<CompressedSA<TText, TSpec, TConfig> const, Rooted>::Type(compressedSA, length(compressedSA));
}

template <typename TText, typename TSpec, typename TConfig>
inline typename Iterator<CompressedSA<TText, TSpec, TConfig>, Standard>::Type
end(CompressedSA<TText, TSpec, TConfig> & compressedSA, Standard const & /* dummy */)
{
    return typename Iterator<CompressedSA<TText, TSpec, TConfig>, Standard>::Type(compressedSA, length(compressedSA));
}

template <typename TText, typename TSpec, typename TConfig>
inline typename Iterator<CompressedSA<TText, TSpec, TConfig> const, Standard>::Type
end(CompressedSA<TText, TSpec, TConfig> const & compressedSA, Standard/* dummy */)
{
    return typename Iterator<CompressedSA<TText, TSpec, TConfig> const, Standard>::Type(compressedSA, length(compressedSA));
}

}
#endif // INDEX_FM_COMPRESSED_SA_ITERATOR_H_
