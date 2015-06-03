// ==========================================================================
//                               misc_memset.h
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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef INCLUDE_SEQAN_MISC_MEMSET_H_
#define INCLUDE_SEQAN_MISC_MEMSET_H_

#include <seqan/basic/basic_metaprogramming.h>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function memset()
// ----------------------------------------------------------------------------

/*!
 * @fn memset
 * @headerfile <seqa/misc/memset.h>
 * @brief An implementation of <tt>memset</tt> with fixed number of bytes using Metaprogramming.
 *
 * @signature template <unsigned SIZE>
 *            void memset(ptr[, c]);
 *
 * @tparam SIZE The number of bytes to reset, i.e. the values in <tt>[ptr, ptr + SIZE)</tt> are set to <tt>c</tt>.
 *
 * @param[in,out] ptr The first position to reset.
 * @param[in]     c   The value to set to.
 */

// TODO(holtgrew): Does memset() really belong in this header? Used in find_myers_ukknonen.h, pump_lcp_core.h, pipe_sample.h, file_async

using std::memset;

// Implementation of memset() with fill size.

template <unsigned SIZE, bool direct>
struct MemsetWorker
{
    finline static
    void run(unsigned char * ptr, unsigned char c)
    {
        std::memset(ptr, c, SIZE);
    }
};

template <unsigned  SIZE>
struct MemsetWorker<SIZE, true>
{
    finline static
    void run(unsigned char* ptr, unsigned char c)
    {
        *((unsigned*)ptr) = ((unsigned)c << 24) + ((unsigned)c << 16) + ((unsigned)c << 8) + (unsigned)c;
        MemsetWorker<SIZE - 4, true>::run(ptr + 4, c);
    }
};

template <>
struct MemsetWorker<0, true>
{
    finline static void
    run(unsigned char*, unsigned char)
    {}
};

template <>
struct MemsetWorker<1, true>
{
    finline static
    void run(unsigned char* ptr, unsigned char c)
    {
        *ptr = c;
    }
};

template <>
struct MemsetWorker<2, true>
{
    finline static
    void run(unsigned char* ptr, unsigned char c)
    {
        *(unsigned short *)ptr = ((unsigned short)c << 8) + (unsigned short)c;
    }
};

template <>
struct MemsetWorker<3, true> {
    finline static
    void run(unsigned char* ptr, unsigned char c)
    {
        MemsetWorker<2, true>::run(ptr, c);
        MemsetWorker<1, true>::run(ptr + 2, c);
    }
};

template <unsigned SIZE>
finline void memset(void* ptr, unsigned char c)
{
    MemsetWorker<SIZE, SIZE <= 32>::run((unsigned char*)ptr, c);
}

// Implementation of memset() with fill value.

template <unsigned SIZE, bool direct, unsigned char c>
struct MemsetConstValueWorker
{
    finline static void run(unsigned char* ptr) { std::memset(ptr, c, SIZE); }
};

template <unsigned  SIZE, unsigned char c>
struct MemsetConstValueWorker<SIZE, true, c>
{
    finline static
    void run(unsigned char* ptr)
    {
        *((unsigned*)ptr) = ((unsigned)c << 24) + ((unsigned)c << 16) + ((unsigned)c << 8) + (unsigned)c;
        MemsetConstValueWorker<SIZE - 4, true, c>::run(ptr + 4);
    }
};

template <unsigned char c>
struct MemsetConstValueWorker<0, true, c>
{
    finline static
    void run(unsigned char*) {}
};

template <unsigned char c>
struct MemsetConstValueWorker<1, true, c>
{
    finline static
    void run(unsigned char* ptr) {
        *ptr = c;
    }
};

template <unsigned char c>
struct MemsetConstValueWorker<2, true, c>
{
    finline static
    void run(unsigned char* ptr)
    {
        *(unsigned short *)ptr = ((unsigned short)c << 8) + (unsigned short)c;
    }
};

template <unsigned char c>
struct MemsetConstValueWorker<3, true, c>
{
    finline static
    void run(unsigned char* ptr)
    {
        MemsetConstValueWorker<2, true, c>::run(ptr);
        MemsetConstValueWorker<1, true, c>::run(ptr + 2);
    }
};

template <unsigned SIZE, unsigned char c>
finline void
memset(void* ptr)
{
    MemsetConstValueWorker<SIZE, SIZE <= 32, c>::run((unsigned char*)ptr);
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_MISC_MEMSET_H_
