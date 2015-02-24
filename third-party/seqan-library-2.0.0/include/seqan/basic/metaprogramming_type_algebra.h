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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// Metaprogramming for types algebra.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_TYPE_ALGEBRA_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_TYPE_ALGEBRA_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename T1, typename T2, typename TSpec> struct Pair;
template <typename TTag, typename TSubList> struct TagList;

// ==========================================================================
// Metafunctions
// ==========================================================================

// ----------------------------------------------------------------------------
// Metafunction Sum<T1, T2>
// ----------------------------------------------------------------------------
// T1 + T2

template <typename T1 = void, typename T2 = void>
struct Sum
{
    typedef TagList<T1, TagList<T2, void> >   Type;
};

// ----------------------------------------------------------------------------
// Metafunction Sum<T1, void>
// ----------------------------------------------------------------------------

template <typename T1>
struct Sum<T1, void>
{
    typedef T1 Type;
};

// ----------------------------------------------------------------------------
// Metafunction Sum<void, T2>
// ----------------------------------------------------------------------------

template <typename T2>
struct Sum<void, T2>
{
    typedef T2 Type;
};

// ----------------------------------------------------------------------------
// Metafunction Sum<void, void>
// ----------------------------------------------------------------------------

//template <>
//struct Sum<void, void>
//{
//    typedef void Type;
//};

// ----------------------------------------------------------------------------
// Metafunction Sum<TTag1, TTagList2>
// ----------------------------------------------------------------------------
// T1 + L2[0..0] = T1 + L2[0]

template <typename TTag1, typename TTag2>
struct Sum<TTag1, TagList<TTag2, void> >
{
    typedef typename Sum<TTag1, TTag2>::Type    Type;
};

// ----------------------------------------------------------------------------
// Metafunction Sum<TTag1, TTagList2>
// ----------------------------------------------------------------------------
// T1 + L2[0..m] = T1 + (L2[0] + .. + L2[m])

template <typename TTag1, typename TTag2, typename TSubList2>
struct Sum<TTag1, TagList<TTag2, TSubList2> >
{
    typedef TagList<TTag1, TagList<TTag2, TSubList2> >  Type;
};

// ----------------------------------------------------------------------------
// Metafunction Sum<TTagList1, TTagList2>
// ----------------------------------------------------------------------------
// L1[0..0] + L2[0..m] = L1[0] + (L2[0] + .. + L2[m])

template <typename TTag1, typename TTag2, typename TSubList2>
struct Sum<TagList<TTag1, void>, TagList<TTag2, TSubList2> >
{
    typedef TagList<TTag1, TagList<TTag2, TSubList2> >  Type;
};

// ----------------------------------------------------------------------------
// Metafunction Sum<TTagList1, TTagList2>
// ----------------------------------------------------------------------------
// L1[0..n] + L2[0..m] = L1[0] + (L1[1..n] + (L2[0] + .. + L2[m]))

template <typename TTag1, typename TSubList1, typename TTag2, typename TSubList2>
struct Sum<TagList<TTag1, TSubList1>, TagList<TTag2, TSubList2> >
{
    typedef TagList<TTag1, typename Sum<TSubList1, TagList<TTag2, TSubList2> >::Type> Type;
};

// ----------------------------------------------------------------------------
// Metafunction Product<T1, T2>
// ----------------------------------------------------------------------------
// T1 x T2

template <typename T1 = void, typename T2 = void>
struct Product
{
    typedef Pair<T1, T2, void>  Type;
};

// ----------------------------------------------------------------------------
// Metafunction Product<T1, void>
// ----------------------------------------------------------------------------

//template <typename T1>
//struct Product<T1, void>
//{
//    typedef void    Type;
//};

// ----------------------------------------------------------------------------
// Metafunction Product<void, T2>
// ----------------------------------------------------------------------------

//template <typename T2>
//struct Product<void, T2>
//{
//    typedef void    Type;
//};

// ----------------------------------------------------------------------------
// Metafunction Product<void, void>
// ----------------------------------------------------------------------------

//template <>
//struct Product<void, void>
//{
//    typedef void    Type;
//};

// ----------------------------------------------------------------------------
// Metafunction Product<TTag1, TTagList2>
// ----------------------------------------------------------------------------
// T1 x L2[0..0] = T1 x L2[0]

template <typename TTag1, typename TTag2>
struct Product<TTag1, TagList<TTag2, void> >
{
    typedef TagList<typename Product<TTag1, TTag2>::Type, void> Type;
};

// ----------------------------------------------------------------------------
// Metafunction Product<TTagList1, TTag2>
// ----------------------------------------------------------------------------
// L1[0..0] x T2 = L1[0] x T2

//template <typename TTag1, typename TTag2>
//struct Product<TagList<TTag1, void>, TTag2>
//{
//    typedef TagList<typename Product<TTag1, TTag2>::Type>   Type;
//};

// ----------------------------------------------------------------------------
// Metafunction Product<TTagList1, TTagList2>
// ----------------------------------------------------------------------------
// L1[0..0] x L2[0..0] = L1[0] x L2[0]

//template <typename TTag1, typename TTag2>
//struct Product<TagList<TTag1, void>, TagList<TTag2, void> >
//{
//    typedef TagList<typename Product<TTag1, TTag2>::Type>   Type;
//};

// ----------------------------------------------------------------------------
// Metafunction Product<TTag1, TTagList2>
// ----------------------------------------------------------------------------
// T1 x L2[0..m] = T1 x L2[0] + T1 x L2[1..m]

template <typename TTag1, typename TTag2, typename TSubList2>
struct Product<TTag1, TagList<TTag2, TSubList2> >
{
    typedef TagList<typename Product<TTag1, TTag2>::Type,
                    typename Product<TTag1, TSubList2>::Type>   Type;
};

// ----------------------------------------------------------------------------
// Metafunction Product<TTagList1, TTag2>
// ----------------------------------------------------------------------------
// L1[0..n] x T2 = L1[0] x T2 + L1[1..n] x T2

//template <typename TTag1, typename TSubList1, typename TTag2>
//struct Product<TagList<TTag1, TSubList1>, TTag2>
//{
//    typedef TagList<typename Product<TTag1, TTag2>::Type,
//                    typename Product<TSubList1, TTag2>::Type>   Type;
//};

// ----------------------------------------------------------------------------
// Metafunction Product<TTagList1, TTagList2>
// ----------------------------------------------------------------------------
// L1[0..0] x L2[0..m] = L1[0] x L2[0] + L1[0] x L2[1..m]

template <typename TTag1, typename TTag2, typename TSubList2>
struct Product<TagList<TTag1, void>, TagList<TTag2, TSubList2> >
{
    typedef TagList<typename Product<TTag1, TTag2>::Type,
                    typename Product<TTag1, TSubList2>::Type> Type;
};

// ----------------------------------------------------------------------------
// Metafunction Product<TTagList1, TTagList2>
// ----------------------------------------------------------------------------
// L1[0..n] x L2[0..0] = L1[0] x L2[0] + L1[1..n] x L2[0]

//template <typename TTag1, typename TSubList1, typename TTag2>
//struct Product<TagList<TTag1, TSubList1>, TagList<TTag2, void> >
//{
//    typedef TagList<typename Product<TTag1, TTag2>::Type,
//                    typename Product<TSubList1, TTag2>::Type> Type;
//};

// ----------------------------------------------------------------------------
// Metafunction Product<TTagList1, TTagList2>
// ----------------------------------------------------------------------------
// L1[0..n] x L2[0..m] = L1[0] x L2[0..m] + L1[1..n] x L2[0..m]

template <typename TTag1, typename TSubList1, typename TTag2, typename TSubList2>
struct Product<TagList<TTag1, TSubList1>, TagList<TTag2, TSubList2> >
{
    typedef typename Sum<typename Product<TTag1, TagList<TTag2, TSubList2> >::Type,
                         typename Product<TSubList1, TagList<TTag2, TSubList2> >::Type>::Type Type;
};

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_TYPE_ALGEBRA_H_
