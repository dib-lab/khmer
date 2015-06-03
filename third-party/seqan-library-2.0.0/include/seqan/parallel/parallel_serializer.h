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
// Class to aquire items, process them in parallel and serialize them in the
// order of acquirement.
// ==========================================================================

#ifndef SEQAN_PARALLEL_PARALLEL_SERIALIZER_H_
#define SEQAN_PARALLEL_PARALLEL_SERIALIZER_H_

namespace seqan {

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Struct SerializerItem
// ----------------------------------------------------------------------------

template <typename TValue>
struct SerializerItem
{
    TValue          val;
    SerializerItem  *next;
    bool            ready;
};

// ----------------------------------------------------------------------------
// Class Serializer
// ----------------------------------------------------------------------------

template <typename TValue, typename TWorker>
class Serializer
{
public:
    typedef SerializerItem<TValue>          TItem;
    typedef TItem *                         TItemPtr;
    typedef ResourcePool<TItem>             TPool;
    typedef typename Size<Serializer>::Type TSize;

    CriticalSection     cs;
    TWorker             worker;
    TItemPtr            first;
    TItemPtr            last;
    TPool               pool;
    bool                stop;

    Serializer() :
        first(NULL),
        last(NULL),
        stop(false)
    {}

    template <typename TArg>
    explicit
    Serializer(TArg &arg, TSize maxItems = 1024) :
        worker(arg),
        first(NULL),
        last(NULL),
        pool(maxItems),
        stop(false)
    {}

    ~Serializer()
    {
        while (first != NULL)
        {
            TItemPtr item = first;
            first = first->next;
            delete item;
        }
    }

    operator bool()
    {
        return !stop;
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TValue, typename TWorker>
inline void
clear(Serializer<TValue, TWorker> & me)
{
    me.stop = false;
    while (me.first != NULL)
    {
        TValue *item = me.first;
        me.first = me.first->next;
        releaseValue(me.recycled, item);
    }
    me.last = NULL;
}

// ----------------------------------------------------------------------------
// Function aquireValue()
// ----------------------------------------------------------------------------

// this function is not thread-safe as it would make
// not much sense to order a stream by the random
// order of executition behind a mutex
template <typename TValue, typename TWorker>
inline TValue *
aquireValue(Serializer<TValue, TWorker> & me)
{
    typedef SerializerItem<TValue> TItem;

    TItem *item = aquireValue(me.pool);
    item->next = NULL;
    item->ready = false;

    // add item to the end of our linked list
    {
        ScopedLock<CriticalSection> lock(me.cs);
        if (me.first == NULL)
            me.first = item;
        else
            me.last->next = item;
        me.last = item;
    }
    return &item->val;
}

// ----------------------------------------------------------------------------
// Function releaseValue()
// ----------------------------------------------------------------------------

template <typename TValue, typename TWorker>
inline bool
releaseValue(Serializer<TValue, TWorker> & me, TValue *ptr)
{
    typedef SerializerItem<TValue> TItem;

    TItem *item = reinterpret_cast<TItem *>(ptr);
    SEQAN_ASSERT_NOT(item->ready);

    // changing me.first or the ready flag must be done synchronized (me.mutex)
    // the thread who changed me.first->ready to be true has to write it.

    // change our ready flag and test if me.first->ready became true
    {
        ScopedLock<CriticalSection> lock(me.cs);
        item->ready = true;
        if (item != me.first)
            return true;
    }

    // ok, if we are here it seems that we are responsible for writing the buffer

    SEQAN_ASSERT(me.first != NULL);

    bool success;
    do
    {
        // process item
        success = me.worker(item->val);

        // remove item from linked list
        {
            ScopedLock<CriticalSection> lock(me.cs);
            me.first = item->next;

            // recycle released items
            releaseValue(me.pool, item);

            // can we leave?
            item = me.first;
            if (item == NULL || !item->ready)
                return success;
        }

        // we continue to process the next buffer
    }
    while (success);

    return false;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_PARALLEL_PARALLEL_SERIALIZER_H_
