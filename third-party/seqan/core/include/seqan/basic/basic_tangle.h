// ==========================================================================
//                               basic_tangle.h
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

#ifndef CORE_INCLUDE_SEQAN_BASIC_BASIC_TANGLE_H_
#define CORE_INCLUDE_SEQAN_BASIC_BASIC_TANGLE_H_

namespace seqan {

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

///.Metafunction.Key.param.T.type:Class.Pair

template <typename TKey, typename TObject, typename TSpec>
struct Key<Pair<TKey, TObject, TSpec> > 
{
    typedef TKey Type;
};

// ----------------------------------------------------------------------------
// Metafunction Cargo
// ----------------------------------------------------------------------------

// TODO(holtgrew): Is this part of some adaption?

///.Metafunction.Cargo.param.T.type:Class.Pair

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

/**
.Function.assignQualities
..cat:Alphabets
..summary:Assign quality values between strings.
..signature:assignQualities(target, source)
..param.target:Target string
...type:nolink:@Class.String@ of any alphabet with qualities, e.g. @Spec.DnaQ@, @Spec.Dna5Q@
..param.source:Source string.
...type:nolink:@Class.String@ of $int$ or $char$.
..remarks:
The target is resized to the length of source.
This function calls @Function.assignQualityValue@ for all entries of $target$ and $source$, look at the documentation of @Function.assignQualityValue@ on how the values of $source$ are interpreted.
..remarks:
Note that qualities are expected to be in PHRED scale.
..see:Function.assignQualityValue
..include:seqan/basic.h
*/

// TODO(holtgrew): Uncomment, place somewhere that knows both iterators and assignQualityValue.
template <typename TDest, typename TSource>
void assignQualities(TDest &dst, TSource const &src)
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

template <typename T>
inline T
unknownValueImpl(T *)
{
    SEQAN_CHECKPOINT;
    return 'N';
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BASIC_BASIC_TANGLE_H_
