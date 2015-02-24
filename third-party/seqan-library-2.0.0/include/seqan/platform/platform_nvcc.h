// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
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
//     * Neither the name of NVIDIA Corporation nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

/*!
 * @macro PLATFORM_CUDA
 * @headerfile <seqan/platform.h>
 * @brief Defined if the compiler is nvcc (or CUDA-capable).
 *
 * @signature #define PLATFORM_CUDA
 *
 * @section Remarks
 *
 * This macro is a synonym for __CUDACC__.
 */

/*!
 * @macro SEQAN_HOST_DEVICE
 * @headerfile <seqan/platform.h>
 * @brief Prefix for functions working both on host and device side.
 *
 * @signature #define SEQAN_HOST_DEVICE
 *
 * This macro can be placed in front of CUDA-compatible functions that can be callable both on host and device side.
 * The macro expands to <tt>__host__ __device__</tt> on CUDA-capable compilers and is ignored otherwise.
 *
 * @see SEQAN_HOST
 * @see SEQAN_DEVICE
 * @see SEQAN_GLOBAL
 *
 * @section Example
 *
 * @code{.cpp}
 * inline SEQAN_HOST_DEVICE void foo(int & x)
 * {
 *     // I can run on the CPU and on the GPU, yay!
 *     x = 10;
 * }
 * @endcode
 */

/*!
 * @macro SEQAN_HOST
 * @headerfile <seqan/platform.h>
 * @brief Prefix for functions working only on host side.
 *
 * @signature #define SEQAN_HOST
 *
 * This macro can be placed in front of functions that can be callable only from host side.
 * The macro expands to <tt>__host__</tt> on CUDA-capable compilers and is ignored otherwise.
 *
 * @see SEQAN_HOST_DEVICE
 * @see SEQAN_DEVICE
 * @see SEQAN_GLOBAL
 */

/*!
 * @macro SEQAN_DEVICE
 * @headerfile <seqan/platform.h>
 * @brief Prefix for functions working only on device side.
 *
 * @signature #define SEQAN_DEVICE
 *
 * This macro must be placed in front of functions that can be callable only from device side.
 * The macro expands to <tt>__device__</tt> on CUDA-capable compilers and is ignored otherwise.
 *
 * @see SEQAN_HOST_DEVICE
 * @see SEQAN_HOST
 * @see SEQAN_GLOBAL
 *
 * @section Remarks
 *
 * Note that a device function containing CUDA intrinsics will not compile on non CUDA-capable compilers. Therefore, to
 * insure graceful compilation, it is still necessary to wrap CUDA-intrinsic code inside __CUDA_ARCH__ defines.
 */

/*!
 * @macro SEQAN_GLOBAL
 * @headerfile <seqan/platform.h>
 * @brief Prefix for functions defining entry points for device side kernels.
 *
 * @signature #define SEQAN_GLOBAL
 *
 * This macro must be placed in front of functions defining entry points for device side kernels.
 * The macro expands to <tt>__global__</tt> on CUDA-capable compilers and is ignored otherwise.
 *
 * @see SEQAN_HOST_DEVICE
 * @see SEQAN_HOST
 * @see SEQAN_DEVICE
 */

#ifndef PLATFORM_CUDA_H_
#define PLATFORM_CUDA_H_

#ifdef __CUDACC__

#define PLATFORM_CUDA

#define SEQAN_FUNC inline __host__ __device__
#define SEQAN_HOST_DEVICE __host__ __device__
#define SEQAN_HOST __host__
#define SEQAN_DEVICE __device__
#define SEQAN_GLOBAL __global__

// NOTE(esiragusa): this solves a problem with nvcc using gcc4.x on MacOS X
#if defined(PLATFORM_GCC) && defined (__APPLE__)
#undef _GLIBCXX_USE_INT128
#endif

#else

#define SEQAN_FUNC inline
#define SEQAN_HOST_DEVICE
#define SEQAN_HOST
#define SEQAN_DEVICE
#define SEQAN_GLOBAL

#endif  // #ifdef __CUDACC__

#endif  // #ifndef PLATFORM_CUDA_H_
