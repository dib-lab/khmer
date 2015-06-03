// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Module header for the sequence module.
//
// The sequence module contains:
//  * Class String and some specializations.
//  * Class StringSet and its specialization.
//  * Adaptions of STL classes to the SeqAn Sequence concept.
//  * Misc sequence-related code such as lexical comparison helpers.
// ==========================================================================

#ifndef SEQAN_HEADER_SEQUENCE_H
#define SEQAN_HEADER_SEQUENCE_H

// ===========================================================================
// Prerequisites.
// ===========================================================================

#include <seqan/basic.h>
#include <seqan/misc/bit_twiddling_functors.h>
#include <seqan/misc/bit_twiddling.h>

// ----------------------------------------------------------------------------
// STL prerequisites.
// ----------------------------------------------------------------------------
// The classes std::string, std::list and std::vector are adapted in this module.

#include <cassert>

#include <map>  // used in string set
#include <string>
#include <list>
#include <vector>
#include <algorithm>
#include <numeric>

// ----------------------------------------------------------------------------
// Thrust prerequisites.
// ----------------------------------------------------------------------------

#ifdef PLATFORM_CUDA
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#endif

// ===========================================================================
// Forwards.
// ===========================================================================

#include <seqan/sequence/sequence_forwards.h>

// ===========================================================================
// Miscellaneous sequence-related code.
// ===========================================================================

#include <seqan/sequence/sequence_lexical.h>

// ===========================================================================
// Segments: Suffixes, Infixes, Prefixes.
// ===========================================================================

#include <seqan/sequence/sequence_interface.h>
#include <seqan/sequence/segment_base.h>
#include <seqan/sequence/segment_infix.h>
#include <seqan/sequence/segment_suffix.h>
#include <seqan/sequence/segment_prefix.h>
#include <seqan/sequence/segment_utils.h>

// ===========================================================================
// Strings.
// ===========================================================================

#include <seqan/sequence/string_base.h>
#include <seqan/sequence/string_array.h>
#include <seqan/sequence/string_alloc.h>
#include <seqan/sequence/string_cstyle.h>
#include <seqan/sequence/string_block.h>
#include <seqan/sequence/string_packed.h>

#include <seqan/sequence/sequence_shortcuts.h>

// ===========================================================================
// Lightweight Containers.
// ===========================================================================

#include <seqan/sequence/container_view.h>

// ===========================================================================
// Adaptions.
// ===========================================================================

#include <seqan/sequence/adapt_std_list.h>
#include <seqan/sequence/adapt_std_string.h>
#include <seqan/sequence/adapt_std_vector.h>
#include <seqan/sequence/adapt_array_pointer.h>
#ifdef PLATFORM_CUDA
#include <seqan/sequence/adapt_thrust_vector.h>
#endif

// ===========================================================================
// StringSets.
// ===========================================================================

#include <seqan/sequence/iter_concat_virtual.h>
#include <seqan/sequence/sequence_concatenator.h>
#include <seqan/sequence/string_set_base.h>
#include <seqan/sequence/string_set_view.h>
#ifdef PLATFORM_CUDA
#include <seqan/sequence/string_set_device.h>
#endif
#include <seqan/sequence/string_set_concat_direct.h>
#include <seqan/sequence/string_set_dependent_tight.h>
#include <seqan/sequence/string_set_dependent_generous.h>
#include <seqan/sequence/string_set_owner.h>
#include <seqan/sequence/string_set_segment.h>

#endif //#ifndef SEQAN_HEADER_...
