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
// Author: Stephan Aiche <stephan.aiche@fu-berlin.de>
// ==========================================================================
// Facade header for module arg_parse.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ARG_PARSE_H_
#define SEQAN_INCLUDE_SEQAN_ARG_PARSE_H_

// ===========================================================================
// Prerequisites.
// ===========================================================================

#include <seqan/basic.h>
//#include <seqan/file.h>
#include <seqan/stream.h>
#include <seqan/sequence.h>

#include <sys/types.h>
#include <sys/stat.h>
#ifndef PLATFORM_WINDOWS
#include <unistd.h>
#endif

// ===========================================================================
// The ArgParse Headers.
// ===========================================================================

#include <seqan/arg_parse/arg_parse_type_support.h>
#include <seqan/arg_parse/xml_support.h>
#include <seqan/arg_parse/tool_doc.h>
#include <seqan/arg_parse/arg_parse_exceptions.h>

#include <seqan/arg_parse/arg_parse_argument.h>
#include <seqan/arg_parse/arg_parse_option.h>

#include <seqan/arg_parse/argument_parser.h>

#include <seqan/arg_parse/arg_parse_doc.h>
#include <seqan/arg_parse/arg_parse_ctd_support.h>
#include <seqan/arg_parse/arg_parse_parse.h>


#endif  // SEQAN_INCLUDE_SEQAN_ARG_PARSE_H_
