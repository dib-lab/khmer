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
// Facade header for basic_iterator submodule.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_BASIC_ITERATOR_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_BASIC_ITERATOR_H_

// --------------------------------------------------------------------------
// Dependencies
// --------------------------------------------------------------------------

#include <seqan/platform.h>
#include <seqan/basic/basic_fundamental.h>
#include <seqan/basic/basic_concept.h>
#include <seqan/basic/basic_alphabet.h>

#include <seqan/basic/basic_type.h>  // TODO(holtgrew): Temporary dependency

// --------------------------------------------------------------------------
// Sub Module Headers
// --------------------------------------------------------------------------

// Iterator interface definition.
#include <seqan/basic/iterator_interface.h>

// Iter base class.
#include <seqan/basic/iterator_base.h>

// Iterator concept.
#include <seqan/basic/iterator_concept.h>

// PropertyMap concept.
#include <seqan/basic/property_map_concept.h>

// Container concept.
#include <seqan/basic/container_concept.h>

// Counting iterator.
#include <seqan/basic/iterator_counting.h>

// Positional iterator.
#include <seqan/basic/iterator_position.h>

// Makes a container out of begin/end iterators
#include <seqan/basic/iterator_range.h>

// Iterator for adapting iterators to Rooted Iterators.
#include <seqan/basic/iterator_adaptor.h>

// Adaption between STL iterators and SeqAn iterators.
#include <seqan/basic/iterator_adapt_std.h>

// TODO(holtgrew): Does proxy really belong here?
#include <seqan/basic/proxy_base.h>
#include <seqan/basic/proxy_iterator.h>

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_BASIC_ITERATOR_H_
