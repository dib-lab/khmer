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
// Utility macros for parallelism.
// ==========================================================================

#ifndef SEQAN_PARALLEL_PARALLEL_MACROS_H_
#define SEQAN_PARALLEL_PARALLEL_MACROS_H_

/**
.Macro.SEQAN_OMP_PRAGMA
..summary:Portable conditional $#pragma$ issuing if OpenMP is enabled.
..cat:Parallelism
..signature:SEQAN_OMP_PRAGMA(x)
..param.x:The string to issue behind $#pragma omp$.
..remarks:This macro uses portable pragma generation, dependent on the macro $_OPENMP$ being defined (as by the OpenMP standard).
..remarks:This is useful for disabling OpenMP pragmas on compilers that do not support OpenMP to suppress warnings.
..example.text:Parallelize loop with OpenMP if OpenMP is enabled:
..example.code:
SEQAN_OMP_PRAGMA(parallel for) // becomes: #pragma omp parallel for
for (int i = 0; i < x; ++i) {
    // Do work.
}
..example.text:Make an addition atomic if OpenMP is enabled:
..example.code:
SEQAN_OMP_PRAGMA(parallel atomic) // becomes: #pragma omp parallel atomic
i += 1;
 */

#ifdef _OPENMP
  #include <omp.h>
  #if defined(PLATFORM_WINDOWS_MINGW) || defined(PLATFORM_GCC)
    // GCC _Pragma operator
    #define SEQAN_DO_PRAGMA(x) _Pragma(#x)
    #define SEQAN_OMP_PRAGMA(x) SEQAN_DO_PRAGMA(omp x)
  #else  // #if defined(PLATFORM_WINDOWS_MINGW) || defined(PLATFORM_GCC)
    // MSVC __pragma-operator
    #define SEQAN_OMP_PRAGMA(x) __pragma (omp x)
  #endif // #if defined(PLATFORM_WINDOWS_MINGW) || defined(PLATFORM_GCC)
#else  // #ifdef _OPENMP
  #define SEQAN_OMP_PRAGMA(x)

  // low-level OpenMP runtime compatibility
  inline void omp_set_num_threads(int) {}
  inline int  omp_get_num_threads()    { return 1; }
  inline int  omp_get_max_threads()    { return 1; }
  inline int  omp_get_thread_num()     { return 0; }

#endif  // #ifdef _OPENMP

#endif  // SEQAN_PARALLEL_PARALLEL_MACROS_H_
