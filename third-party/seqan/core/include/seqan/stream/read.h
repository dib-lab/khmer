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
// Main File for Record and Document-Reading. Contains only doc right now.
// ==========================================================================


#ifndef SEQAN_STREAM_READ_H_
#define SEQAN_STREAM_READ_H_

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
.Function.readRecord
..class:Class.RecordReader
..cat:Input/Output
..summary:reads one record (e.g. a single DNA-sequence and its meta data) from a @Concept.StreamConcept@, by the means of @Class.RecordReader@
..signature:readRecord(<format specific>, TRecordReader & reader, TTag const &)
..param.<format specific>: possibly multiple fields (e.g. meta and sequence)
..param.reader:The reader object to read from
...type:Class.RecordReader
..param.TTag:The file format tag
..remarks: If not noted otherwise, only a Single-Pass implementation is available for a the given format
..see:Function.read2
..include:seqan/stream.h
*/

/**
.Function.read2
..class:Class.RecordReader
..cat:Input/Output
..summary:reads an entire document from a @Concept.StreamConcept@, by the means of @Class.RecordReader@
..signature:read2(<format specific>, TRecordReader & reader, TTag const &)
..param.<format specific>: possibly multiple StringSets (e.g. of meta and sequences)
..param.reader:The reader object to read from
...type:Class.RecordReader
..param.TTag:The file format tag
..status:Should be renamed to "read" once the old IO-Code is removed
..remarks: This is only supported for Double-Pass IO. If you cannot use Double-Pass IO (e.g. when you cannot seek on the stream), loop over @Function.readRecord@ instead.
..remarks:If not noted otherwise an especially efficient version of the function is used if all StringSets are specialized as @Spec.ConcatDirect@ -StringSets.
..include:seqan/stream.h
..see:Class.RecordReader
..see:Function.readRecord
*/

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_READ_FASTA_FASTQ_H_
