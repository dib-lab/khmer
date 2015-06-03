// ==========================================================================
//                               basic_tangle.h
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

// TODO(holtgrew): This header contains code that does not clearly belongs somewhere else.

#ifndef INCLUDE_SEQAN_BASIC_BASIC_TANGLE_H_
#define INCLUDE_SEQAN_BASIC_BASIC_TANGLE_H_

namespace seqan {

// TODO(holtgrew): Remove this define.
#define SEQAN_NAMESPACE_MAIN seqan

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// TODO(holtgrew): Remove auto-sequence!
template <typename TValue, typename TSpec>
struct Iterator<SimpleType<TValue, TSpec>, Standard>
{
    typedef SimpleType<TValue, TSpec> * Type;
//  typedef Iter<SimpleType<TValue, TSpec>, SimpleIterator> * Type;
};

template <typename TValue, typename TSpec>
struct Iterator<SimpleType<TValue, TSpec> const, Standard>
{
    typedef SimpleType<TValue, TSpec> const * Type;
//  typedef Iter<SimpleType<TValue, TSpec> const, SimpleIterator> * Type;
};

// ----------------------------------------------------------------------------
// Metafunction Key
// ----------------------------------------------------------------------------

// TODO(holtgrew): Is this part of some adaption?

template <typename TKey, typename TObject, typename TSpec>
struct Key<Pair<TKey, TObject, TSpec> >
{
    typedef TKey Type;
};

// ----------------------------------------------------------------------------
// Metafunction Cargo
// ----------------------------------------------------------------------------

// TODO(holtgrew): Is this part of some adaption?

template <typename TKey, typename TCargo, typename TSpec>
struct Cargo<Pair<TKey, TCargo, TSpec> >
{
    typedef TCargo Type;
};


// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function assignQualities()
// ----------------------------------------------------------------------------

/*!
 * @fn assignQualities
 * @headerfile <seqan/basic.h>
 * @brief Assign quality values between strings.
 *
 * @signature void assignQualities(target, source);
 *
 * @param[out] target Target string, can be a String of DnaQ or Dna5Q characters.
 * @param[in]  source Source string.  Can be a String of int or char.
 *
 * @section Remarks
 *
 * The target is resized to the length of source.  This function calls assignQualityValue for all entries of
 * <tt>target</tt> and <tt>source</tt>, look at the documentation of assignQualityValue on how the values of
 * <tt>source</tt> are interpreted.
 *
 * Note that qualities are expected to be in PHRED scale.
 *
 * @see AlphabetWithQualitiesConcept#assignQualityValue
 */

template <typename TDest, typename TSource>
inline void
_assignQualities(TDest &dst, TSource const &src, True)
{
    typedef typename Iterator<TDest>::Type TDestIter;
    typedef typename Iterator<TSource const>::Type TSourceIter;

    if (length(dst) < length(src))
        resize(dst, length(src));

    TDestIter itDst = begin(dst, Standard());
    TSourceIter itSrcEnd = end(src, Standard());

    for (TSourceIter itSrc = begin(src, Standard()); itSrc != itSrcEnd; ++itDst, ++itSrc)
        assignQualityValue(*itDst, *itSrc);
}

template <typename TDest, typename TSource>
inline void
_assignQualities(TDest &, TSource const &, False)
{}

// TODO(holtgrew): Uncomment, place somewhere that knows both iterators and assignQualityValue, maybe in module sequence?
template <typename TDest, typename TSource>
inline void
assignQualities(TDest &dst, TSource const &src)
{
    typedef typename Value<TDest>::Type TValue;
    _assignQualities(dst, src, typename Or<IsSameType<TValue, DnaQ>, IsSameType<TValue, Dna5Q> >::Type());
}

template <typename T>
inline T
unknownValueImpl(T *)
{
    SEQAN_CHECKPOINT;
    return 'N';
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_BASIC_BASIC_TANGLE_H_
