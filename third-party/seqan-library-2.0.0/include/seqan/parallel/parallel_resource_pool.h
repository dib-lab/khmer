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
// Pool for reusable items, e.g. buffer strings
// ==========================================================================

#ifndef SEQAN_PARALLEL_PARALLEL_RESOURCE_POOL_H_
#define SEQAN_PARALLEL_PARALLEL_RESOURCE_POOL_H_

namespace seqan {

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class ResourcePool
// ----------------------------------------------------------------------------

template <typename TValue>
struct ResourcePool
{
    typedef ConcurrentQueue<TValue *, Suspendable<> >   TStack;
    typedef typename Size<TStack>::Type                 TSize;

    TStack recycled;

    ResourcePool(TSize maxSize)
    {
        setWriterCount(recycled, 1);
        for (; maxSize != 0; --maxSize)
            appendValue(recycled, (TValue *)NULL);
    }

    ~ResourcePool()
    {
        unlockWriting(recycled);
        TValue *ptr = NULL;
        unsigned count = 0;
        while (popBack(ptr, recycled))
        {
            if (ptr != NULL)
                count++;
            delete ptr;
        }
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function aquireValue()
// ----------------------------------------------------------------------------

template <typename TValue>
inline TValue *
aquireValue(ResourcePool<TValue> & me)
{
    TValue *ptr = NULL;
    if (!popBack(ptr, me.recycled))
        return NULL;

    if (ptr == NULL)
        ptr = new TValue;

    return ptr;
}

// ----------------------------------------------------------------------------
// Function releaseValue()
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
releaseValue(ResourcePool<TValue> & me, TValue *ptr)
{
    appendValue(me.recycled, ptr);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_PARALLEL_PARALLEL_RESOURCE_POOL_H_
