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
// Facade header for module bam_io.
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_H_

// ===========================================================================
// Prerequisites.
// ===========================================================================

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/align.h>
#include <seqan/misc/misc_name_store_cache.h>

// ===========================================================================
// Data Structures & Conversion.
// ===========================================================================

#include <seqan/bam_io/bam_io_context.h>
#include <seqan/bam_io/cigar.h>
#include <seqan/bam_io/bam_alignment_record.h>
#include <seqan/bam_io/bam_header_record.h>
#include <seqan/bam_io/bam_sam_conversion.h>
#include <seqan/bam_io/bam_tags_dict.h>

// ===========================================================================
// Actual I/O Code.
// ===========================================================================

#include <seqan/bam_io/read_sam.h>
#include <seqan/bam_io/write_sam.h>

// BAM I/O is only available when ZLIB is available (and thus the BGz Stream).
#if SEQAN_HAS_ZLIB
#include <seqan/bam_io/read_bam.h>
#include <seqan/bam_io/write_bam.h>
#endif  // #if SEQAN_HAS_ZLIB

// ===========================================================================
// BAM Index Related.
// ===========================================================================

// BAM indices are only available when ZLIB is available.
#if SEQAN_HAS_ZLIB
#include <seqan/bam_io/bam_index_base.h>
#include <seqan/bam_io/bam_index_bai.h>
#endif  // #if SEQAN_HAS_ZLIB

// ===========================================================================
// Utility Routines.
// ===========================================================================

#include <seqan/bam_io/bam_alignment_record_util.h>

// ===========================================================================
// Easy BAM / SAM I/O.
// ===========================================================================

#include <seqan/bam_io/xam_reader.h>
#if SEQAN_HAS_ZLIB
#include <seqan/bam_io/bam_reader.h>
#endif  // #if SEQAN_HAS_ZLIB
#include <seqan/bam_io/sam_reader.h>

#include <seqan/bam_io/xam_writer.h>
#if SEQAN_HAS_ZLIB
#include <seqan/bam_io/bam_writer.h>
#endif  // #if SEQAN_HAS_ZLIB
#include <seqan/bam_io/sam_writer.h>

#include <seqan/bam_io/bam_stream.h>

#endif  // CORE_INCLUDE_SEQAN_BAM_IO_H_
