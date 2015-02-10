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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// Main File for record- and Document-writing. Contains only doc right now.
// ==========================================================================


#ifndef SEQAN_STREAM_WRITE_H_
#define SEQAN_STREAM_WRITE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

/**
.Function.writeRecord
..cat:Input/Output
..summary:write one record (e.g. a single DNA-sequence and its meta data) to a @Concept.StreamConcept@
..signature:writeRecord(TStream & stream, <format specific>, TTag const &)
..param.stream:The Stream object to write to
...type:Concept.StreamConcept
..param.format specific: possibly multiple fields (e.g. meta and sequence)
..param.TTag:The file format tag
..see:Function.write2
..include:seqan/stream.h
*/

/**
.Function.write2
..cat:Input/Output
..summary:writes an entire document to a @Concept.StreamConcept@
..signature:write2(TStream & stream, <format specific>, TTag const &)
..param.stream:The Stream object to write to
...type:Concept.StreamConcept
..param.format specific: possibly multiple StringSets (e.g. of meta and sequences)
..param.TTag:The file format tag
..status:Should be renamed to "write" once the old IO-Code is removed
..include:seqan/stream.h
..see:Function.writeRecord
*/

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_WRITE_H_
