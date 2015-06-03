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

/*!
 * @macro PLATFORM_WINDOWS
 * @headerfile <seqan/platform.h>
 * @brief Defined if the compiler targets Windows (MSVC or MinGW).
 *
 * @signature #define PLATFORM_WINDOWS
 *
 * @macro PLATFORM_WINDOWS_MINGW
 * @headerfile <seqan/platform.h>
 * @brief Defined if the compiler is MinGW.
 *
 * @signature #define PLATFORM_WINDOWS_MINGW
 */

#ifndef PLATFORM_WINDOWS
#define PLATFORM_WINDOWS

#define PLATFORM_WINDOWS_MINGW

#define finline __inline__

// The symbols SEQAN_IS_64_BIT and SEQAN_IS_32_BIT can be used to check
// whether we are on a 32 bit or on a 64 bit machine.
#if defined(__amd64__) || defined(__x86_64__) || defined(__ia64__)
#define SEQAN_IS_64_BIT 1
#define SEQAN_IS_32_BIT 0
#else  // #if defined(__amd64__) || defined(__x86_64__) || defined(__ia64__)
#define SEQAN_IS_64_BIT 0
#define SEQAN_IS_32_BIT 1
#endif  // #if defined(__amd64__) || defined(__x86_64__) || defined(__ia64__)

#include <inttypes.h>

// Define unsigned variants of builtin Windows compiler types.
typedef unsigned __int64 __uint64;
typedef unsigned __int32 __uint32;
typedef unsigned __int16 __uint16;
typedef unsigned __int8 __uint8;

// Define ftello
#ifndef ftello
#define ftello(fp) ftell(fp)
#endif  // #ifndef ftello

// detect C++11 support
#if defined(__GXX_EXPERIMENTAL_CXX0X__)
#  define SEQAN_CXX11_STANDARD
#endif

#endif  // #ifndef PLATFORM_WINDOWS
