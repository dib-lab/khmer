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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Facade header for module seq_io.
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_SEQ_IO_H_
#define CORE_INCLUDE_SEQAN_SEQ_IO_H_

// ===========================================================================
// Prerequisites.
// ===========================================================================

#include <seqan/basic.h>
#include <seqan/stream.h>

// ===========================================================================
// Lower-Level I/O Interface for Sequences
// ===========================================================================

// File Format Guessing
#include <seqan/seq_io/guess_stream_format.h>

#include <seqan/seq_io/read_fasta_fastq.h>
#include <seqan/seq_io/read_embl.h>
#include <seqan/seq_io/read_genbank.h>

#include <seqan/seq_io/write_fasta_fastq.h>

// ===========================================================================
// SequenceStream
// ===========================================================================

#include <seqan/seq_io/sequence_stream_impl.h>
#include <seqan/seq_io/sequence_stream.h>

// ===========================================================================
// Reading FASTA for Demos
// ===========================================================================

#include <seqan/seq_io/simple_read_fasta.h>

// ===========================================================================
// Genomic Region
// ===========================================================================

#include <seqan/seq_io/genomic_region.h>

// ===========================================================================
// FAI Index
// ===========================================================================

#include <seqan/seq_io/fai_index.h>

#endif  // CORE_INCLUDE_SEQAN_SEQ_IO_H_
