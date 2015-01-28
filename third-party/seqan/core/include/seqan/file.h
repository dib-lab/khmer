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

#ifndef SEQAN_HEADER_FILE_H
#define SEQAN_HEADER_FILE_H

//____________________________________________________________________________
// prerequisites

#include <iostream>
#include <climits>
#include <cstdio>
#include <list>
#include <vector>
#include <map>
#include <cmath>

#include <seqan/sequence.h>
#include <seqan/modifier.h>

#ifdef SEQAN_NEW_IO
#include <cctype>
#endif


//____________________________________________________________________________

#include <seqan/file/file_forwards.h>

#ifndef SEQAN_NEW_IO
#include <seqan/file/cstream.h>
#include <seqan/file/stream.h>
#endif

#include <seqan/file/chunk_collector.h>
// #include <seqan/file/meta.h>

//____________________________________________________________________________
// files

#include <seqan/file/file_interface.h>
#ifndef SEQAN_NEW_IO
#include <seqan/file/file_cstyle.h>
#endif
#include <seqan/file/file_base.h>

#include <seqan/system.h>	// async file (default file type of File<>)

/*#include <seqan/system/file_sync.h>
#include <seqan/system/system_event.h>
#include <seqan/system/file_async.h>
*/

//____________________________________________________________________________
// file formats
#ifndef SEQAN_NEW_IO
#include <seqan/file/file_filereaderiterator.h>
#include <seqan/file/file_filereader.h>

#include <seqan/file/file_format.h>

#include <seqan/file/stream_algorithms.h>

//file formats for sequences
#include <seqan/file/file_format_raw.h>
#include <seqan/file/file_format_fasta.h>
#include <seqan/file/file_format_embl.h>
#include <seqan/file/file_format_genbank.h>

//file formats for alignments
#include <seqan/file/file_format_fasta_align.h>

//others
#include <seqan/file/file_format_cgviz.h>
#endif

//____________________________________________________________________________

//#include <seqan/file/file_format_guess.h>

//____________________________________________________________________________
// external strings

#include <seqan/file/file_page.h>
#include <seqan/file/file_mapping.h>
#include <seqan/file/string_mmap.h>
//#include <seqan/file/file_pager.h>
#include <seqan/file/string_external.h>

#ifndef SEQAN_NEW_IO
#include <seqan/file/file_format_mmap.h>
#endif

#endif //#ifndef SEQAN_HEADER_...
