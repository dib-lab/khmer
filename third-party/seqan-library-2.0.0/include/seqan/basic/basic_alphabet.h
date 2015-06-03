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
// Facade header for the basic/alphabet sub module.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_BASIC_ALPHABET_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_BASIC_ALPHABET_H_

// --------------------------------------------------------------------------
// Dependencies
// --------------------------------------------------------------------------

#include <seqan/platform.h>
#include <seqan/basic/basic_metaprogramming.h>
#include <seqan/basic/basic_fundamental.h>
#include <seqan/basic/basic_concept.h>

// --------------------------------------------------------------------------
// Sub Module Headers
// --------------------------------------------------------------------------

// The first, more generic part of this sub module contains concepts, and the
// necessary metafunction and function prototypes as well as some forward and
// default implementations.

// The alphabet concepts.
#include <seqan/basic/alphabet_concept.h>

// Alphabet math functions / metafunctions forwards and default
// implementations.
#include <seqan/basic/alphabet_math.h>

// Adaptions for builtin C++ types to the alphabet concepts.
// TODO(holtgrew): Move into second part?
#include <seqan/basic/alphabet_adapt_builtins.h>

// Alphabets from Bioinformatics application - gapped unknown values.
#include <seqan/basic/alphabet_bio.h>

// Alphabet with qualities metafunction and function forwards and default
// implementations
#include <seqan/basic/alphabet_qualities.h>

// Forward declarations and default implementations for storage related
// aspects of alphabets.
#include <seqan/basic/alphabet_storage.h>

// The second part contains basic alphabets, namely the SimpleType class with
// the specialization for the biological data types, such as Dna, Rna, and
// Amino Acid.  Furthermore, it contains a character storing profile entries.

// The SimpleType class.
#include <seqan/basic/alphabet_simple_type.h>

// Conversion tables for the biological SimpleType specializations.
#include <seqan/basic/alphabet_residue_tabs.h>

// Conversion functions for the biological SimpleType specializations.
#include <seqan/basic/alphabet_residue_funcs.h>

// The actual biological SimpleType specializations.
#include <seqan/basic/alphabet_residue.h>

// The profile character implementation.
#include <seqan/basic/alphabet_profile.h>

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_BASIC_ALPHABET_H_
