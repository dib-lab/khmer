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

#ifndef SEQAN_HEADER_FILE_H
#define SEQAN_HEADER_FILE_H

// ===========================================================================
// Prerequisites.
// ===========================================================================

// ----------------------------------------------------------------------------
// STL prerequisites
// ----------------------------------------------------------------------------

#include <iostream>
#include <climits>
#include <cstdio>
#include <list>
#include <vector>
#include <map>
#include <cmath>
#include <sstream>
#include <cctype>

// ----------------------------------------------------------------------------
// SeqAn prerequisites
// ----------------------------------------------------------------------------

#include <seqan/sequence.h>
#include <seqan/modifier.h>

// ===========================================================================
// Forwards.
// ===========================================================================

#include <seqan/file/file_forwards.h>

// ===========================================================================
// File.
// ===========================================================================

#include <seqan/file/file_interface.h>
#include <seqan/file/file_base.h>
#include <seqan/system.h>    // async file (default file type of File<>)

// ===========================================================================
// External Strings.
// ===========================================================================

// ----------------------------------------------------------------------------
// Paging and mapping
// ----------------------------------------------------------------------------

#include <seqan/file/file_page.h>
#include <seqan/file/file_mapping.h>
//#include <seqan/file/file_pager.h>

// ----------------------------------------------------------------------------
// String specializations
// ----------------------------------------------------------------------------

#include <seqan/file/string_mmap.h>
#include <seqan/file/string_external.h>

#endif //#ifndef SEQAN_HEADER_...
