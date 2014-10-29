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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_SYSTEM_BASE_H
#define SEQAN_HEADER_SYSTEM_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

#if SEQAN_ENABLE_DEBUG  // Note the new-style for macros, is always defined and 0/1

#define SEQAN_DO_SYS(_cond) SEQAN_ASSERT(_cond)
#define SEQAN_DO_SYS1(_cond) SEQAN_DO_SYS(_cond)
#define SEQAN_DO_SYS2(_cond, _comment) SEQAN_ASSERT_MSG(_cond, _comment)

#else  // #ifdef SEQAN_ENABLE_DEBUG

#if defined(PLATFORM_GCC) || defined(PLATFORM_WINDOWS_MINGW)
// GCC warns below that the "value computed is not used".  However,
// MSVC does not like casting void values to void. Thus, this
// distinction.
#define SEQAN_DO_SYS(_cond) do { (void) _cond; } while (false)
#else   // #if defined(PLATFORM_GCC) || defined(PLATFORM_WINDOWS_MINGW)
#define SEQAN_DO_SYS(_cond) do { _cond; } while (false)
#endif  // #if defined(PLATFORM_GCC) || defined(PLATFORM_WINDOWS_MINGW)

#define SEQAN_DO_SYS1(_cond) SEQAN_DO_SYS(_cond)
#define SEQAN_DO_SYS2(_cond, _comment) SEQAN_DO_SYS(_cond)

#endif  // #ifdef SEQAN_ENABLE_DEBUG

}

#endif
