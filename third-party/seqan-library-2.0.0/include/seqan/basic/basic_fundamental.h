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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Facade hader for the submodule basic_fundamental.
//
// This module contains fundamental code such as forward declarations and
// prototypes for common metafunctions like Value<>, functions assign() etc.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_BASIC_FUNDAMENTAL_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_BASIC_FUNDAMENTAL_H_

// --------------------------------------------------------------------------
// Prerequisites
// --------------------------------------------------------------------------

#include <algorithm>

#include <seqan/platform.h>
#include <seqan/basic/basic_metaprogramming.h>

// --------------------------------------------------------------------------
// Sub Module Headers
// --------------------------------------------------------------------------

// Macros for deprecating code.
#include <seqan/basic/macro_deprecated.h>

// Pseudo header with documentation for builtin functions.
#include <seqan/basic/builtin_functions.h>

// Common metafunctions such as Value<>.
#include <seqan/basic/fundamental_metafunctions.h>

// Basic tag-related code.
#include <seqan/basic/fundamental_tags.h>

// Functions and metafunctions to use contiguous chunks of memory
#include <seqan/basic/fundamental_chunking.h>

// Definition of assign(), set(), move().
#include <seqan/basic/fundamental_transport.h>

// Code supporting comparison.
#include <seqan/basic/fundamental_comparison.h>

// Conversion support.
#include <seqan/basic/fundamental_conversion.h>

// TODO(holtgrew): This is not fundamental.  Should go into sequence module.
// Construct/destruct functions for arrays.
#include <seqan/basic/array_construct_destruct.h>

// TODO(holtgrew): This is not really fundamental, either. Should go into its own sub module.
// Hosted type.
#include <seqan/basic/hosted_type_interface.h>

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_HOSTED_TYPE_INTERFACE_H_
