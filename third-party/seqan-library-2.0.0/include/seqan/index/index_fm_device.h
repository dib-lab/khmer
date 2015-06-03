// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
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
//     * Neither the name of NVIDIA Corporation nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE
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
// This file contains Index specializations for thrust::device_vector.
// ==========================================================================

#ifndef INDEX_FM_DEVICE_H_
#define INDEX_FM_DEVICE_H_

namespace seqan {

// ----------------------------------------------------------------------------
// Typedef DnaStringSet
// ----------------------------------------------------------------------------

typedef StringSet<DnaString, Owner<ConcatDirect<> > >   DnaStringSet;

// ----------------------------------------------------------------------------
// Typedef CudaFMIndexConfig
// ----------------------------------------------------------------------------

struct CudaFMIndexConfig
{
    typedef TwoLevels<void>    TValuesSpec;
    typedef Naive<void>        TSentinelsSpec;

    static const unsigned SAMPLING = 10;
};

typedef FMIndex<void, CudaFMIndexConfig>        CudaFMIndexSpec;

// ----------------------------------------------------------------------------
// Typedef DnaStringSetFMIndex
// ----------------------------------------------------------------------------

typedef Index<DnaString, CudaFMIndexSpec>       DnaStringFMIndex;
typedef Index<DnaStringSet, CudaFMIndexSpec>    DnaStringSetFMIndex;

// ----------------------------------------------------------------------------
// Metafunction SAValue
// ----------------------------------------------------------------------------

template <>
struct SAValue<DnaStringSet>
{
    typedef Pair<__uint8, __uint32, Pack> Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size; Index
// ----------------------------------------------------------------------------

template <>
struct Size<DnaStringSetFMIndex>
{
    typedef __uint32 Type;
};

template <>
struct Size<View<DnaStringSetFMIndex>::Type>
{
    typedef __uint32 Type;
};

template <>
struct Size<Device<DnaStringSetFMIndex>::Type>
{
    typedef __uint32 Type;
};

template <>
struct Size<View<Device<DnaStringSetFMIndex>::Type>::Type>
{
    typedef __uint32 Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size; LF
// ----------------------------------------------------------------------------

template <>
struct Size<LF<DnaStringSet, void, CudaFMIndexConfig> >
{
    typedef __uint32 Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size; Rank Dictionary
// ----------------------------------------------------------------------------

template <typename TSpec>
struct Size<RankDictionary<Dna, TwoLevels<TSpec> > >
{
    typedef __uint32 Type;
};

template <typename TSpec>
struct Size<RankDictionary<bool, TwoLevels<TSpec> > >
{
    typedef __uint32 Type;
};

template <typename TSpec>
struct Size<RankDictionary<bool, Naive<TSpec> > >
{
    typedef __uint32 Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value; Shape
// ----------------------------------------------------------------------------

template <unsigned q>
struct Value<Shape<Dna, UngappedShape<q> > >
{
    typedef __uint32    Type;
};

}

#endif  // #ifndef INDEX_FM_DEVICE_H_
