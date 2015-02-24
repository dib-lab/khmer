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

#ifndef SEQAN_STRING_SET_VIEW_H
#define SEQAN_STRING_SET_VIEW_H

namespace seqan {

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction View
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec>
struct View<StringSet<TString, TSpec> >
{
    typedef StringSet<typename View<TString>::Type, TSpec>      Type;
};

// ----------------------------------------------------------------------------
// Metafunction RemoveView
// ----------------------------------------------------------------------------

template <typename TString, typename TViewSpec, typename TSpec>
struct RemoveView<StringSet<ContainerView<TString, TViewSpec>, TSpec> >
{
    typedef StringSet<TString, TSpec>   Type;
};

// ----------------------------------------------------------------------------
// Metafunction IsView
// ----------------------------------------------------------------------------

template <typename TString, typename TViewSpec, typename TSpec>
struct IsView<StringSet<ContainerView<TString, TViewSpec>, TSpec> > : public True {};

// ----------------------------------------------------------------------------
// Metafunction StringSetLimits
// ----------------------------------------------------------------------------

template <typename TString, typename TViewSpec, typename TSpec>
struct StringSetLimits<StringSet<ContainerView<TString, TViewSpec>, TSpec> >
{
    typedef typename View<typename StringSetLimits<StringSet<TString, TSpec> >::Type>::Type     Type;
};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function _initStringSetLimits()
// --------------------------------------------------------------------------

template <typename TString, typename TViewSpec, typename TSpec>
inline void
_initStringSetLimits(StringSet<ContainerView<TString, TViewSpec>, TSpec> & /* me */) {}


}  // namespace seqan

#endif  // #ifndef SEQAN_STRING_SET_VIEW_H
