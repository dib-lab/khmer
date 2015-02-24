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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// 2-level spinlock and corresponding scoped locks for each level
// ==========================================================================

#ifndef SEQAN_PARALLEL_PARALLEL_LOCK_H_
#define SEQAN_PARALLEL_PARALLEL_LOCK_H_

#if defined(__SSE2__) && !defined(__CUDACC__)
#include <xmmintrin.h>  // _mm_pause()
#endif

#ifdef PLATFORM_WINDOWS
#include <Windows.h>
#else
#include <sched.h>
#endif

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

struct Mutex;
inline void yieldProcessor();

// ============================================================================
// Classes
// ============================================================================

class SpinDelay
{
public:
    static const unsigned LOOPS_BEFORE_YIELD = 16;
    unsigned duration;

    SpinDelay() :
        duration(1)
    {}
};

inline void
clear(SpinDelay & me)
{
    me.duration = 1;
}

inline void
waitFor(SpinDelay & me)
{
    if (me.duration <= me.LOOPS_BEFORE_YIELD)
    {
        for (unsigned i = me.duration; i != 0; --i)
            yieldProcessor();
        me.duration *= 2;
    }
    else
    {
#ifdef PLATFORM_WINDOWS
#if _WIN32_WINNT >= 0x0400
        SwitchToThread();
#endif
#else
        sched_yield();
#endif
    }
}

template <typename TAtomic, typename TValue>
inline void
spinWhileEq(TAtomic & x, TValue cmp)
{
    SpinDelay spinDelay;
    while (x == cmp)
        waitFor(spinDelay);
}

template <typename TAtomic, typename TValue>
inline void
spinWhileNeq(TAtomic & x, TValue cmp)
{
    SpinDelay spinDelay;
    while (x != cmp)
        waitFor(spinDelay);
}

template <typename TAtomic, typename TValue>
inline void
spinCas(TAtomic & x, TValue cmp, TValue y)
{
    SpinDelay spinDelay;
#ifdef SEQAN_CXX11_STL
    TValue exp = cmp;
    while (!x.compare_exchange_weak(exp, y))
    {
        exp = cmp;
        waitFor(spinDelay);
    }
#else
    while (!atomicCasBool(x, cmp, y))
        waitFor(spinDelay);
#endif
}

// ----------------------------------------------------------------------------
// Class ReadWriteLock
// ----------------------------------------------------------------------------
// this lock augments a class by thread-safety as follows:
//  - supports multiple concurrent readers (possibly waiting for writer to finish)
//  - supports only a single writer at a time (possibly waiting for readers or other writers to finish)
//  - the writer has higher priority than all readers

class ReadWriteLock
{
public:
    Atomic<unsigned>::Type readers;
    Atomic<unsigned>::Type writers;

    ReadWriteLock() :
        readers(0),
        writers(0)
    {}
};

// ----------------------------------------------------------------------------
// Class ScopedLock
// ----------------------------------------------------------------------------

template <typename TMutex = Mutex, typename TParallel = Parallel>
struct ScopedLock
{
    TMutex & mutex;

    explicit
    ScopedLock(TMutex & mutex) :
        mutex(mutex)
    {
        lock(mutex);
    }

    ~ScopedLock()
    {
        unlock(mutex);
    }

};

template <typename TLock>
struct ScopedLock<TLock, Serial>
{
    explicit
    ScopedLock(TLock &)
    {}
};

// ----------------------------------------------------------------------------
// Class ScopedReadLock
// ----------------------------------------------------------------------------

template <typename TLock = ReadWriteLock, typename TParallel = Parallel>
struct ScopedReadLock
{
    TLock & lock;

    explicit
    ScopedReadLock(TLock & lock) :
        lock(lock)
    {
        lockReading(lock);
    }

    ~ScopedReadLock()
    {
        unlockReading(lock);
    }

};

template <typename TLock>
struct ScopedReadLock<TLock, Serial>
{
    explicit
    ScopedReadLock(TLock &)
    {}
};

// ----------------------------------------------------------------------------
// Class ScopedWriteLock
// ----------------------------------------------------------------------------

template <typename TLock = ReadWriteLock, typename TParallel = Parallel>
struct ScopedWriteLock
{
    TLock & lock;

    explicit
    ScopedWriteLock(TLock & lock) :
        lock(lock)
    {
        lockWriting(lock);
    }

    ~ScopedWriteLock()
    {
        unlockWriting(lock);
    }

};

template <typename TLock>
struct ScopedWriteLock<TLock, Serial>
{
    explicit
    ScopedWriteLock(TLock &)
    {}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function yieldProcessor()
// ----------------------------------------------------------------------------

inline void
yieldProcessor()
{
#if defined( __CUDACC__)
    // don't wait on the GPU
#elif defined(PLATFORM_WINDOWS_VS)
    YieldProcessor();
#elif defined(__SSE2__)
    _mm_pause();
#else
    __asm__ __volatile__("rep; nop" : : );
#endif
}

// ----------------------------------------------------------------------------
// Function lockReading()
// ----------------------------------------------------------------------------

inline void
lockReading(ReadWriteLock & lock)
{
    do
    {
        // wait for the end of a write access
        spinWhileNeq(lock.writers, 0u);

        atomicInc(lock.readers);

        if (lock.writers == 0u)
            break;

        // writer hasn't noticed us -> retry
        atomicDec(lock.readers);
    }
    while (true);
}

// ----------------------------------------------------------------------------
// Function unlockReading()
// ----------------------------------------------------------------------------

inline void
unlockReading(ReadWriteLock & lock)
{
    atomicDec(lock.readers);
}

// ----------------------------------------------------------------------------
// Function lockWriting()
// ----------------------------------------------------------------------------

inline void
lockWriting(ReadWriteLock & lock)
{
    // wait until we are the only writer
    spinCas(lock.writers, 0u, 1u);

    // wait until all readers are done
    spinWhileNeq(lock.readers, 0u);
}

// ----------------------------------------------------------------------------
// Function unlockWriting()
// ----------------------------------------------------------------------------

inline void
unlockWriting(ReadWriteLock & lock)
{
    lock.writers = 0;
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

inline bool
empty(ReadWriteLock & lock)
{
    return (lock.readers == 0u && lock.writers == 0u);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_PARALLEL_PARALLEL_LOCK_H_
