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

#ifndef PLATFORM_GCC
#define PLATFORM_GCC

// should be set before including anything
#ifndef _FILE_OFFSET_BITS
  #define _FILE_OFFSET_BITS 64
#endif

#ifndef _LARGEFILE_SOURCE
  #define _LARGEFILE_SOURCE
#endif

// The symbols SEQAN_IS_64_BIT and SEQAN_IS_32_BIT can be used to check
// whether we are on a 32 bit or on a 64 bit machine.
#if defined(__amd64__) || defined(__x86_64__) || defined(__ia64__)
#define SEQAN_IS_64_BIT 1
#define SEQAN_IS_32_BIT 0
#else  // #if defined(__amd64__) || defined(__x86_64__) || defined(__ia64__)
#define SEQAN_IS_64_BIT 0
#define SEQAN_IS_32_BIT 1
#endif  // #if defined(__amd64__) || defined(__x86_64__) || defined(__ia64__)

//#include <unistd.h>
#include <inttypes.h>

#define finline __inline__

// default 64bit type
typedef uint64_t __uint64; // nolint

// default 32bit type
typedef int32_t __int32;   // nolint
typedef uint32_t __uint32; // nolint

// default 16bit type
typedef int16_t __int16;   // nolint
typedef uint16_t __uint16; // nolint

// default 8bit type
typedef int8_t __int8;     // nolint
typedef uint8_t __uint8;   // nolint

// detect gcc C++11 support
#if defined(__GXX_EXPERIMENTAL_CXX0X__)
#  define SEQAN_CXX11_STANDARD
#endif

// detect clang C++11 support
#ifdef __has_feature
#  if __has_feature(cxx_static_assert)
#    define SEQAN_CXX11_STANDARD
#  endif
#endif

#endif  // #ifndef PLATFORM_GCC
