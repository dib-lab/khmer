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

#ifndef SEQAN_STRING_SET_DEVICE_H
#define SEQAN_STRING_SET_DEVICE_H

namespace seqan {

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Device                                              [StringSet]
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec>
struct Device<StringSet<TString, TSpec> >
{
    typedef StringSet<typename Device<TString>::Type, TSpec>    Type;
};

// ----------------------------------------------------------------------------
// Metafunction IsDevice                                     [Device StringSet]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlloc, typename TSpec>
struct IsDevice<StringSet<thrust::device_vector<TValue, TAlloc>, TSpec> > : public True {};

// ----------------------------------------------------------------------------
// Metafunction StringSetLimits                              [Device StringSet]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlloc, typename TSpec>
struct StringSetLimits<StringSet<thrust::device_vector<TValue, TAlloc>, TSpec> >
{
    typedef thrust::device_vector<TValue, TAlloc>   TString_;
    typedef typename Size<TString_>::Type           TSize_;
    typedef thrust::device_vector<TSize_>           Type;
};


}  // namespace seqan

#endif  // #ifndef SEQAN_STRING_SET_DEVICE_H
