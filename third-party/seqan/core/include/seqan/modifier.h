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
// Umbrella header for the modifier module.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_MODIFIER_H_
#define SEQAN_CORE_INCLUDE_SEQAN_MODIFIER_H_

// --------------------------------------------------------------------------
// Prerequisites
// --------------------------------------------------------------------------

#include <functional>
#include <list>

#include <seqan/basic.h>
#include <seqan/sequence.h>

// --------------------------------------------------------------------------
// Modified Alphabets
// --------------------------------------------------------------------------

#include <seqan/modifier/modifier_alphabet.h>
#include <seqan/modifier/modifier_alphabet_expansion.h>

// --------------------------------------------------------------------------
// Modified Strings and Iterators Base Classes
// --------------------------------------------------------------------------

#include <seqan/modifier/modifier_iterator.h>
#include <seqan/modifier/modifier_string.h>

// --------------------------------------------------------------------------
// Applications of Base Classes
// --------------------------------------------------------------------------

#include <seqan/modifier/modifier_functors.h>
#include <seqan/modifier/modifier_view.h>
#include <seqan/modifier/modifier_reverse.h>
#include <seqan/modifier/modifier_shortcuts.h>

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_MODIFIER_H_
