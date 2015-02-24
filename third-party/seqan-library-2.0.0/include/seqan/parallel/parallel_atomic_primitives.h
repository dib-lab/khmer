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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

// SEQAN_NO_GENERATED_FORWARDS: No forwards are generated for this file.

#ifndef SEQAN_PARALLEL_PARALLEL_ATOMIC_PRIMITIVES_H_
#define SEQAN_PARALLEL_PARALLEL_ATOMIC_PRIMITIVES_H_

#if defined(PLATFORM_WINDOWS) && !defined(PLATFORM_WINDOWS_MINGW)
#include <intrin.h>
#endif  // #if defined(PLATFORM_WINDOWS) && !defined(PLATFORM_WINDOWS_MINGW)

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

#ifdef SEQAN_CXX11_STL

template <typename T>
struct Atomic
{
    typedef std::atomic<T> Type;
};

#else

template <typename T>
struct Atomic
{
    typedef T Type;
};

#endif

// ============================================================================
// Functions
// ============================================================================

/*!
 * @defgroup AtomicPrimitives Atomic Primitives
 * @brief Portable atomic operations.
 */

/*!
 * @fn AtomicPrimitives#atomicInc
 * @headerfile <seqan/parallel.h>
 * @brief Atomically increment an integer.
 *
 * @signature TResult atomicInc(x);
 *
 * @param[in,out] x An integer, by reference.
 *
 * @return TResult The old value of $x$, <tt>TResult</tt> has the same type as <tt>x</tt>.
 *
 * @section Remarks
 *
 * This is equivalent to an atomic <tt>++x</tt>.
 *
 * Note that atomic increments are limited to 32 bit and 64 bit with MSVC (64 bit is only available on 64 bit Windows).
 *
 * You are responsible for correctly aligning <tt>x</tt> such that the atomic increment works on the hardware you
 * target.
 */

/*!
 * @fn AtomicPrimitives#atomicDec
 * @headerfile <seqan/parallel.h>
 * @brief Atomically decrement an integer.
 *
 * @signature TResult atomicDec(x);
 *
 * @param[in,out] x An integer, by reference.
 *
 * @return TResult The old value of $x$, <tt>TResult</tt> has the same type as <tt>x</tt>.
 *
 * @section Remarks
 *
 * This is equivalent to an atomic <tt>--x</tt>.
 *
 * Note that atomic decrements are limited to 32 bit and 64 bit with MSVC (64 bit is only available on 64 bit Windows).
 *
 * You are responsible for correctly aligning <tt>x</tt> such that the atomic decrement works on the hardware you
 * target.
 */

/*!
 * @fn AtomicPrimitives#atomicAdd
 * @headerfile <seqan/parallel.h>
 * @brief Atomically add an integer to another integer.
 *
 * @signature TResult atomicAdd(x, y)
 *
 * @param[in,out] x Integer, by reference.
 * @param[in]     y Integer to add to <tt>x</tt>.
 *
 * @return TResult The old value of <tt>x</tt>.
 *
 * @section Remarks
 *
 * This is equivalent to an atomic <tt>x += y</tt>.
 *
 * Note that atomic fetch-and-add is limited to 32 bit and 64 bit with MSVC (64 bit is only available on 64 bit
 * Windows).
 *
 * You are responsible for correctly aligning <tt>x</tt> such that the atomic increment works on the hardware you
 * target.
 */

/*!
 * @fn AtomicPrimitives#atomicOr
 * @headerfile <seqan/parallel.h>
 * @brief Atomically combine two integers with <tt>OR</tt> operation.
 *
 * @signature TResult atomicOr(x, y);
 *
 * @param[in,out] x Integer, by reference.
 * @param[in]     y Integer to combine with <tt>OR</tt> operation.
 *
 * @return TResult The old value of <tt>x</tt>, <tt>TResult</tt> is the type of <tt>x</tt>.
 *
 * @section Remarks
 *
 * This is equivalent to an atomic <tt>x |= y</tt>.
 *
 * Atomic fetch-and-or for 64 bit integers is only available on 64 bit processors when targeting Intel.
 *
 * Atomic fetch-and-or does not work in VS8 on 64 bit Windows, you can only use <tt>atomicOr()</tt> portably on 32 and
 * 64 bit integers.
 *
 * You are responsible for correctly aligning <tt>x</tt> such that the atomic increment works on the hardware you
 * target.
 */

/*!
 * @fn AtomicPrimitives#atomicXor
 * @headerfile <seqan/parallel.h>
 * @brief Atomically combine two integers with <tt>XOR</tt> operation.
 *
 * @signature TResult atomicXor(x, y);
 *
 * @param[in,out] x Integer, by reference.
 * @param[in]     y Integer to combine with <tt>XOR</tt> operation.
 *
 * @return TResult The old value of <tt>x</tt>, <tt>TResult</tt> is the type of <tt>x</tt>.
 *
 * @section Remarks
 *
 * This is equivalent to an atomic <tt>x ^= y</tt>.
 *
 * Atomic fetch-and-xor fxor 64 bit integers is only available on 64 bit processxors when targeting Intel.
 *
 * Atomic fetch-and-xor does not wxork in VS8 on 64 bit Windows, you can only use <tt>atomicXor()</tt> pxortably on 32 and
 * 64 bit integers.
 *
 * You are responsible fxor cxorrectly aligning <tt>x</tt> such that the atomic increment wxorks on the hardware you
 * target.
 */

/*!
 * @fn AtomicPrimitives#atomicCas
 * @headerfile <seqan/parallel.h>
 * @brief Atomic ompare-and-Swap operation.
 *
 * @signature TResult atomicCas(x, cmp, y)
 *
 * @param[in,out] x   Pointer to the integer to swap.
 * @param[in,out] cmp Value to compare <tt>x</tt> with.
 * @param[in]     y   Value to set <tt>x</tt> to if it is equal to <tt>cmp</tt>.
 *
 * @return TResult Returns the original value of x.
 *
 * @section Remarks
 *
 * The pseudo code for this is as follows:
 *
 * @code{.cpp}
 * atomic {
 *     T val = *(&x);
 *     if (val == cmp)
 *         *(&x) = y;
 *     return val;
 * }
 * @endcode
 *
 * On Windows, atomic CAS is only available for 16, 32, and 64 bit integers, 64 bit is only available on 64 bit Windows.
 *
 * You are responsible for correctly aligning <tt>x</tt> such that the atomic increment works on the hardware you
 * target.
 */

// TODO(holtgrew): What about correct alignment?!

#ifndef SEQAN_CACHE_LINE_SIZE
#define SEQAN_CACHE_LINE_SIZE 128
#endif

#if defined(PLATFORM_WINDOWS) && !defined(PLATFORM_WINDOWS_MINGW)

// ----------------------------------------------------------------------------
// Implementation in MSVC
// ----------------------------------------------------------------------------

// We break the standard code layout here since we only wrap compiler
// intrinsics and it's easier to see things with one glance this way.
#pragma intrinsic(_InterlockedOr, _InterlockedXor, _InterlockedCompareExchange)

template <typename T, typename S>
inline T _atomicOr(T volatile &x, ConstInt<sizeof(char)>, S y) { return _InterlockedOr8(reinterpret_cast<char volatile *>(&x), y); }
template <typename T, typename S>
inline T _atomicXor(T volatile &x, ConstInt<sizeof(char)>, S y) { return _InterlockedXor8(reinterpret_cast<char volatile *>(&x), y); }

template <typename T>
inline T _atomicInc(T volatile &x, ConstInt<sizeof(short)>) { return InterlockedIncrement16(reinterpret_cast<short volatile *>(&x)); }
template <typename T>
inline T _atomicDec(T volatile &x, ConstInt<sizeof(short)>) { return InterlockedDecrement16(reinterpret_cast<short volatile *>(&x)); }
template <typename T, typename S>
inline T _atomicAdd(T volatile &x, ConstInt<sizeof(short)>, S y) { return InterlockedExchangeAdd16(reinterpret_cast<short volatile *>(&x), y); }
template <typename T, typename S>
inline T _atomicOr(T volatile &x, ConstInt<sizeof(short)>, S y) { return _InterlockedOr16(reinterpret_cast<short volatile *>(&x), y); }
template <typename T, typename S>
inline T _atomicXor(T volatile &x, ConstInt<sizeof(short)>, S y) { return _InterlockedXor16(reinterpret_cast<short volatile *>(&x), y); }
template <typename T, typename S, typename U>
inline T _atomicCas(T volatile &x, ConstInt<sizeof(short)>, S cmp, U y) { return _InterlockedCompareExchange16(reinterpret_cast<short volatile *>(&x), y, cmp); }

template <typename T>
inline T _atomicInc(T volatile &x, ConstInt<sizeof(LONG)>) { return InterlockedIncrement(reinterpret_cast<LONG volatile *>(&x)); }
template <typename T>
inline T _atomicDec(T volatile &x, ConstInt<sizeof(LONG)>) { return InterlockedDecrement(reinterpret_cast<LONG volatile *>(&x)); }
template <typename T, typename S>
inline T _atomicAdd(T volatile &x, ConstInt<sizeof(LONG)>, S y) { return InterlockedExchangeAdd(reinterpret_cast<LONG volatile *>(&x), y); }
template <typename T, typename S>
inline T _atomicOr(T volatile &x, ConstInt<sizeof(long)>, S y) { return _InterlockedOr(reinterpret_cast<long volatile *>(&x), y); }
template <typename T, typename S>
inline T _atomicXor(T volatile &x, ConstInt<sizeof(long)>, S y) { return _InterlockedXor(reinterpret_cast<long volatile *>(&x), y); }
template <typename T, typename S, typename U>
inline T _atomicCas(T volatile &x, ConstInt<sizeof(long)>, S cmp, U y) { return _InterlockedCompareExchange(reinterpret_cast<long volatile *>(&x), y, cmp); }

#ifdef _WIN64
template <typename T>
inline T _atomicInc(T volatile &x, ConstInt<sizeof(LONGLONG)>) { return InterlockedIncrement64(reinterpret_cast<LONGLONG volatile *>(&x)); }
template <typename T>
inline T _atomicDec(T volatile &x, ConstInt<sizeof(LONGLONG)>) { return InterlockedDecrement64(reinterpret_cast<LONGLONG volatile *>(&x)); }
template <typename T, typename S>
inline T _atomicAdd(T volatile &x, ConstInt<sizeof(LONGLONG)>, S y) { return InterlockedExchangeAdd64(reinterpret_cast<LONGLONG volatile *>(&x), y); }
template <typename T, typename S>
inline T _atomicOr(T volatile &x, ConstInt<sizeof(__int64)>, S y) { return _InterlockedOr64(reinterpret_cast<__int64 volatile *>(&x), y); }
template <typename T, typename S>
inline T _atomicXor(T volatile &x, ConstInt<sizeof(__int64)>, S y) { return _InterlockedXor64(reinterpret_cast<__int64 volatile *>(&x), y); }
template <typename T, typename S, typename U>
inline T _atomicCas(T volatile &x, ConstInt<sizeof(__int64)>, S cmp, U y) { return _InterlockedCompareExchange64(reinterpret_cast<__int64 volatile *>(&x), y, cmp); }
#endif  // #ifdef _WIN64

template <typename T>
inline T atomicInc(T volatile & x) { return _atomicInc(x, ConstInt<sizeof(T)>()); }
template <typename T>
inline T atomicDec(T volatile & x) { return _atomicDec(x, ConstInt<sizeof(T)>()); }
template <typename T, typename S>
inline T atomicAdd(T volatile &x, S y) { return _atomicAdd(x, ConstInt<sizeof(T)>(), y); }
template <typename T, typename S>
inline T atomicOr(T volatile &x, S y) { return _atomicOr(x, ConstInt<sizeof(T)>(), y); }
template <typename T, typename S>
inline T atomicXor(T volatile &x, S y) { return _atomicXor(x, ConstInt<sizeof(T)>(), y); }
template <typename T, typename S, typename U>
inline T atomicCas(T volatile &x, S cmp, U y) { return _atomicCas(x, ConstInt<sizeof(T)>(), cmp, y); }
template <typename T, typename S, typename U>
inline bool atomicCasBool(T volatile &x, S cmp, U y) { return _atomicCas(x, ConstInt<sizeof(T)>(), cmp, y) == cmp; }

template <typename T>
inline T atomicPostInc(T volatile & x) { return atomicInc(x) - 1; }
template <typename T>
inline T atomicPostDec(T volatile & x) { return atomicDec(x) + 1; }


#else  // #if defined(PLATFORM_WINDOWS) && !defined(PLATFORM_WINDOWS_MINGW)

// ----------------------------------------------------------------------------
// Implementation in GCC (LLVM is GCC compatible)
// ----------------------------------------------------------------------------

template <typename T>
inline T atomicInc(T volatile & x)
{
    return __sync_add_and_fetch(&x, 1);
}

template <typename T>
inline T atomicPostInc(T volatile & x)
{
    return __sync_fetch_and_add(&x, 1);
}

template <typename T>
inline T atomicDec(T volatile & x)
{
    return __sync_add_and_fetch(&x, -1);
}

template <typename T>
inline T atomicPostDec(T volatile & x)
{
    return __sync_fetch_and_add(&x, -1);
}

template <typename T1, typename T2>
inline T1 atomicAdd(T1 volatile & x, T2 y)
{
    return __sync_add_and_fetch(&x, y);
}

template <typename T>
inline T atomicOr(T volatile & x, T y)
{
    return __sync_or_and_fetch(&x, y);
}

template <typename T>
inline T atomicXor(T volatile & x, T y)
{
    return __sync_xor_and_fetch(&x, y);
}

template <typename T>
inline T atomicCas(T volatile & x, T cmp, T y)
{
    return __sync_val_compare_and_swap(&x, cmp, y);
}

template <typename T>
inline bool atomicCasBool(T volatile & x, T cmp, T y)
{
    return __sync_bool_compare_and_swap(&x, cmp, y);
}

template <typename T>
inline T atomicSwap(T volatile & x, T y)
{
    return __sync_lock_test_and_set(x, y);
}

// Pointer versions

template <typename T>
inline T * atomicInc(T * volatile & x)
{
    return (T *) __sync_add_and_fetch((size_t volatile *)&x, sizeof(T));
}

template <typename T>
inline T * atomicPostInc(T * volatile & x)
{
    return (T *) __sync_fetch_and_add((size_t volatile *)&x, sizeof(T));
}

template <typename T>
inline T * atomicDec(T * volatile & x)
{
    return (T *) __sync_add_and_fetch((size_t volatile *)&x, -sizeof(T));
}

template <typename T>
inline T * atomicPostDec(T * volatile & x)
{
    return (T *) __sync_fetch_and_add((size_t volatile *)&x, -sizeof(T));
}

template <typename T1, typename T2>
inline T1 * atomicAdd(T1 * volatile & x, T2 y)
{
    return (T1 *) __sync_add_and_fetch((size_t volatile *)&x, y * sizeof(T2));
}

#endif  // #if defined(PLATFORM_WINDOWS) && !defined(PLATFORM_WINDOWS_MINGW)


// ----------------------------------------------------------------------------
// Wrappers to use faster non-synced functions in serial implementations
// ----------------------------------------------------------------------------

template <typename T>   inline T atomicInc(T          & x,             Serial)      { return ++x;                    }
template <typename T>   inline T atomicPostInc(T      & x,             Serial)      { return x++;                    }
template <typename T>   inline T atomicDec(T          & x,             Serial)      { return --x;                    }
template <typename T>   inline T atomicPostDec(T      & x,             Serial)      { return x--;                    }
template <typename T>   inline T atomicOr (T          & x, T y,        Serial)      { return x |= y;                 }
template <typename T>   inline T atomicXor(T          & x, T y,        Serial)      { return x ^= y;                 }
// In serial mode, there is no other thread changing cmp except us
template <typename T>   inline T atomicCas(T          & x, T cmp, T y, Serial)      { if (x == cmp) x = y; return x; }
//template <typename T>   inline bool atomicCasBool(T   & x, T cmp, T y, Serial)      { if (x == cmp) { x = y; return true; } return false; }
template <typename T>   inline bool atomicCasBool(T volatile & x, T, T y, Serial)   { x = y; return true;            }

template <typename T>   inline T atomicInc(T volatile & x,             Parallel)    { return atomicInc(x);           }
template <typename T>   inline T atomicPostInc(T volatile & x,         Parallel)    { return atomicPostInc(x);       }
template <typename T>   inline T atomicDec(T volatile & x,             Parallel)    { return atomicDec(x);           }
template <typename T>   inline T atomicPostDec(T volatile & x,         Parallel)    { return atomicPostDec(x);       }
template <typename T>   inline T atomicOr (T volatile & x, T y,        Parallel)    { return atomicOr(x, y);         }
template <typename T>   inline T atomicXor(T volatile & x, T y,        Parallel)    { return atomicXor(x, y);        }
template <typename T>   inline T atomicCas(T volatile & x, T cmp, T y, Parallel)    { return atomicCas(x, cmp, y);   }
template <typename T>   inline bool atomicCasBool(T volatile & x, T cmp, T y, Parallel) { return atomicCasBool(x, cmp, y); }

template <typename T1, typename T2>   inline T1 atomicAdd(T1          & x, T2 y, Serial)    { return x = x + y; }
template <typename T1, typename T2>   inline T1 atomicAdd(T1 volatile & x, T2 y, Parallel)  { return atomicAdd(x, y); }


// C++11 atomic wrappers

#ifdef SEQAN_CXX11_STL
template <typename T>   inline T atomicInc(std::atomic<T>        & x     )        { return ++x;                    }
template <typename T>   inline T atomicPostInc(std::atomic<T>    & x     )        { return x++;                    }
template <typename T>   inline T atomicDec(std::atomic<T>        & x     )        { return --x;                    }
template <typename T>   inline T atomicPostDec(std::atomic<T>    & x     )        { return x--;                    }
template <typename T>   inline T atomicOr (std::atomic<T>        & x, T y)        { return x |= y;                 }
template <typename T>   inline T atomicXor(std::atomic<T>        & x, T y)        { return x ^= y;                 }
template <typename T>   inline T atomicCas(std::atomic<T>        & x, T cmp, T y, Serial)   { if (x == cmp) x = y;             return x;   }
template <typename T>   inline T atomicCas(std::atomic<T>        & x, T cmp, T y, Parallel) { x.compare_exchange_weak(cmp, y); return cmp; }
template <typename T>   inline bool atomicCasBool(std::atomic<T> & x, T    , T y, Serial)   { x = y; return true;                          }
template <typename T>   inline bool atomicCasBool(std::atomic<T> & x, T cmp, T y, Parallel) { return x.compare_exchange_weak(cmp, y);      }
#endif  // #ifdef SEQAN_CXX11_STL

} // namespace seqan

#endif  // #if defined(PLATFORM_WINDOWS) && !defined(PLATFORM_WINDOWS_MINGW)
